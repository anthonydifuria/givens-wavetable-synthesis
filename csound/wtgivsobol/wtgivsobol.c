/*
 * wtgivsobol.c
 * ============================================================================
 * Csound plugin opcode: wtgivsobol
 *
 * DESCRIZIONE
 * -----------
 * Oscillatore wavetable con mescolanza timbrica guidata da SEQUENZE DI SOBOL
 * applicate come angoli di rotazioni di Givens.
 *
 * ── Cos'è una sequenza di Sobol? ──────────────────────────────────────────
 * Una sequenza di Sobol è una sequenza quasi-casuale a bassa discrepanza:
 * i punti sono "più uniformemente distribuiti" di punti pseudo-random puri.
 * In 1D, la sequenza Sobol riempie [0,1) in modo molto regolare, evitando
 * cluster e vuoti. Questo produce angoli di rotazione distribuiti in modo
 * quasi-uniforme su [-1,1], garantendo che ogni regione dello spazio timbrico
 * SO(1024) venga esplorata in modo bilanciato.
 *
 * Algoritmo Sobol 1D (Gray-code increment):
 *   x_n = x_{n-1} XOR v_c     dove c = trailing_zeros(n)
 *   v_c = 1 << (31 - c)        (direction numbers semplificati)
 *
 * ── Comportamento DETERMINISTICO ─────────────────────────────────────────
 * IMPORTANTE: il Sobol viene RESETTATO a ogni k-block (seed=0).
 * Questo significa che gli stessi parametri producono SEMPRE la stessa
 * wavetable — comportamento identico alla versione Python.
 * È una mappa stabile parametri → timbro, non una deriva casuale.
 * kmove=0 restituisce esattamente la wavetable base.
 *
 * SINTASSI CSOUND
 * ---------------
 *   aOut wtgivsobol kfreq, kamp, itab,
 *                       kstart, kright, kleft, kboth,
 *                       ktheta0, kdecay, kmove
 *
 * PARAMETRI
 * ---------
 *   kfreq   (k) : frequenza [Hz]
 *   kamp    (k) : ampiezza
 *   itab    (i) : numero ftable (lunghezza >= 1024)
 *   kstart  (k) : indice di partenza [0..1023]
 *   kright  (k) : densità mixing locale forward [0..1]
 *   kleft   (k) : densità mixing locale backward [0..1]
 *   kboth   (k) : densità mixing non-locale (passo primo) [0..1]
 *   ktheta0 (k) : scala angoli [radianti]
 *   kdecay  (k) : decadimento geometrico [0..1]
 *   kmove   (k) : intensità trasformazione [0..1]
 *                 kmove=0 => wavetable base inalterata
 * ============================================================================
 */

#include "csdl.h"
#include <math.h>
#include <string.h>
#include <stdint.h>

#define WT_LEN  1024
#define WT_MASK (WT_LEN - 1)

static inline double clampd(double v, double lo, double hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}
static inline int clampi(int v, int lo, int hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}

/* ============================================================
 * givens_inplace
 * x'[i] =  cos(t)*x[i] - sin(t)*x[j]
 * x'[j] =  sin(t)*x[i] + cos(t)*x[j]
 * ============================================================ */
static inline void givens_inplace(MYFLT *x, int i, int j, double theta) {
  double c = cos(theta), s = sin(theta);
  double xi = (double)x[i], xj = (double)x[j];
  x[i] = (MYFLT)(c*xi - s*xj);
  x[j] = (MYFLT)(s*xi + c*xj);
}

static void peak_normalize(MYFLT *x) {
  double m = 0.0;
  for (int i = 0; i < WT_LEN; i++) {
    double a = fabs((double)x[i]);
    if (a > m) m = a;
  }
  if (m < 1e-12) return;
  double inv = 1.0 / m;
  for (int i = 0; i < WT_LEN; i++) x[i] = (MYFLT)((double)x[i] * inv);
}

/* ============================================================
 * SOBOL 1D — generatore deterministico a bassa discrepanza
 *
 * Stato: contatore intero i + accumulatore x (32 bit XOR).
 * Ad ogni chiamata a sobol_next01():
 *   1) incrementa i
 *   2) calcola c = numero di zeri finali di i (count trailing zeros)
 *   3) direction number: v = 1 << (31 - c)
 *   4) x ^= v
 *   5) ritorna x / 2^32  in [0,1)
 *
 * La sequenza risultante ha discrepanza O(log(N)/N), molto migliore
 * di pseudo-random O(1/sqrt(N)).
 * ============================================================ */
typedef struct {
  uint32_t i;   /* contatore chiamate (parte da 0) */
  uint32_t x;   /* accumulatore XOR (stato corrente) */
} SOBOL1D;

/* Conta i bit zero finali (trailing zeros) di n>0 */
static inline int ctz32(uint32_t n) {
  int c = 0;
  while ((n & 1u) == 0u) { n >>= 1; c++; }
  return c;
}

/* Prossimo campione in [0,1) */
static inline double sobol_next01(SOBOL1D *s) {
  s->i += 1u;
  int c = ctz32(s->i);
  uint32_t v = (uint32_t)(1u << (31 - c));
  s->x ^= v;
  return (double)s->x / 4294967296.0;   /* / 2^32 */
}

/* Reset deterministico: stesso seed => stessa sequenza */
static inline void sobol_reset(SOBOL1D *s) {
  s->i = 0u;
  s->x = 0u;
}

/* ============================================================
 * Struttura dati dell'opcode
 * ============================================================ */
typedef struct {
  OPDS h;
  MYFLT *out;

  MYFLT *kfreq, *kamp, *itab;
  MYFLT *kstart, *kright, *kleft, *kboth, *ktheta0, *kdecay, *kmove;

  double phase;
  double sr;

  MYFLT base[WT_LEN];   /* copia statica della wavetable originale */
  MYFLT wt[WT_LEN];     /* copia di lavoro trasformata */

  SOBOL1D sob;           /* generatore Sobol (resettato ogni k-block) */
} OP;

/* Carica wavetable base dalla ftable Csound */
static int load_base(CSOUND *csound, OP *p) {
  FUNC *ft = csound->FTnp2Find(csound, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) return NOTOK;
  for (int i = 0; i < WT_LEN; i++) p->base[i] = ft->ftable[i];
  return OK;
}

/* ============================================================
 * next_theta
 * Genera il prossimo angolo di rotazione dalla sequenza Sobol.
 *
 * u  ∈ [0,1)  da Sobol
 * z  = 2u-1  ∈ [-1,1]  (centrato su zero)
 * theta = move * theta0 * decay^t * z
 *
 * Il decadimento geometrico decay^t fa sì che le prime rotazioni
 * (vicine a kstart) abbiano angoli più grandi delle successive.
 * ============================================================ */
static inline double next_theta(OP *p, int t, double move, double theta0, double decay) {
  double u = sobol_next01(&p->sob);
  double z = 2.0*u - 1.0;
  return move * theta0 * pow(decay, (double)t) * z;
}

/* ============================================================
 * build_wavetable
 * Costruisce la wavetable trasformata applicando rotazioni Givens
 * con angoli guidati dalla sequenza Sobol.
 *
 * RESET DETERMINISTICO: sobol_reset() viene chiamato all'inizio,
 * quindi la stessa combinazione di parametri produce sempre
 * la stessa wavetable (identico comportamento al codice Python).
 * ============================================================ */
static void build_wavetable(OP *p, int start,
                            double right, double left, double both,
                            double theta0, double decay, double move)
{
  /* Reset Sobol: FONDAMENTALE per comportamento deterministico.
   * Senza reset, il Sobol continuerebbe ad avanzare tra un k-block
   * e l'altro, producendo deriva incontrollata anche a parametri fissi. */
  sobol_reset(&p->sob);

  memcpy(p->wt, p->base, WT_LEN * sizeof(MYFLT));

  const int PRIME_STEP = 37;   /* coprimo con 1024: garantisce visita globale */
  const int MAX_STEPS  = 512;  /* budget CPU per k-block */

  int stepsB = (int)floor(both  * (double)MAX_STEPS);
  int stepsR = (int)floor(right * (double)MAX_STEPS);
  int stepsL = (int)floor(left  * (double)MAX_STEPS);

  int t = 0;   /* contatore globale passi (per il decadimento) */

  /* BOTH: mixing non-locale a passo primo */
  int i = start & WT_MASK;
  for (int s = 0; s < stepsB; s++) {
    int j = (i + PRIME_STEP) & WT_MASK;
    double th = next_theta(p, t++, move, theta0, decay);
    givens_inplace(p->wt, i, j, th);
    i = (i + 1) & WT_MASK;   /* drift per evitare ripetizioni di coppie */
  }

  /* RIGHT: mixing locale forward */
  i = start & WT_MASK;
  for (int s = 0; s < stepsR; s++) {
    int j = (i + 1) & WT_MASK;
    double th = next_theta(p, t++, move, theta0, decay);
    givens_inplace(p->wt, i, j, th);
    i = j;
  }

  /* LEFT: mixing locale backward */
  i = start & WT_MASK;
  for (int s = 0; s < stepsL; s++) {
    int j = (i - 1) & WT_MASK;
    double th = next_theta(p, t++, move, theta0, decay);
    givens_inplace(p->wt, i, j, th);
    i = j;
  }

  peak_normalize(p->wt);
}

/* Loop audio con interpolazione lineare */
static void render(OP *p, MYFLT *out, uint32_t offset, uint32_t nsmps) {
  double sr    = p->sr;
  double freq  = (double)*p->kfreq;  if (freq < 0.0) freq = 0.0;
  double amp   = (double)*p->kamp;
  double phase = p->phase;
  double inc   = (sr > 0.0) ? (freq / sr) : 0.0;

  for (uint32_t n = offset; n < nsmps; n++) {
    phase += inc;
    phase -= floor(phase);

    double idx  = phase * (double)WT_LEN;
    int    i0   = (int)idx;
    double frac = idx - (double)i0;
    int    i1   = (i0 + 1) & WT_MASK;

    out[n] = (MYFLT)(amp * ((double)p->wt[i0] + ((double)p->wt[i1] - (double)p->wt[i0]) * frac));
  }
  p->phase = phase;
}

/* ============================================================
 * op_init  [i-rate]
 * ============================================================ */
static int op_init(CSOUND *csound, OP *p) {
  p->phase = 0.0;
  p->sr    = (double)csound->GetSr(csound);
  sobol_reset(&p->sob);   /* seed deterministico = 0 */
  return OK;
}

/* ============================================================
 * op_perf  [a+k-rate]
 * Ogni k-block:
 *   1) Carica wavetable base
 *   2) Reset Sobol + build_wavetable
 *   3) Render audio
 * ============================================================ */
static int op_perf(CSOUND *csound, OP *p) {
  MYFLT   *out   = p->out;
  uint32_t offset = p->h.insdshead->ksmps_offset;
  uint32_t early  = p->h.insdshead->ksmps_no_end;
  uint32_t nsmps  = CS_KSMPS;

  if (UNLIKELY(offset)) memset(out, 0, offset * sizeof(MYFLT));
  if (UNLIKELY(early))  { nsmps -= early; memset(&out[nsmps], 0, early * sizeof(MYFLT)); }

  if (UNLIKELY(load_base(csound, p) != OK)) {
    for (uint32_t n = offset; n < nsmps; n++) out[n] = FL(0.0);
    return OK;
  }

  int start = clampi((int)floor((double)*p->kstart + 0.5), 0, WT_LEN - 1);

  double right  = clampd((double)*p->kright,  0.0, 1.0);
  double left   = clampd((double)*p->kleft,   0.0, 1.0);
  double both   = clampd((double)*p->kboth,   0.0, 1.0);
  double theta0 = (double)*p->ktheta0;
  double decay  = clampd((double)*p->kdecay,  0.0, 1.0);
  double move   = clampd((double)*p->kmove,   0.0, 1.0);

  build_wavetable(p, start, right, left, both, theta0, decay, move);
  render(p, out, offset, nsmps);
  return OK;
}

/* ============================================================
 * Registrazione opcode
 * Firma: aOut wtgivsobol kfreq, kamp, itab,
 *                            kstart, kright, kleft, kboth,
 *                            ktheta0, kdecay, kmove
 * ============================================================ */
static OENTRY localops[] = {
  { "wtgivsobol", sizeof(OP), 0, 3, "a", "kkikkkkkkk",
    (SUBR)op_init, (SUBR)op_perf }
};

LINKAGE
