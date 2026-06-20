/*
 * wtgivnlocal.c
 * ============================================================================
 * Csound plugin opcode: wtgivnlocal
 *
 * DESCRIZIONE
 * -----------
 * Oscillatore wavetable con mescolanza timbrica NON-LOCALE tramite ROTAZIONI
 * DI GIVENS a passo primo.
 *
 * A differenza di wtgivlocal (che opera su coppie adiacenti), questo opcode
 * usa un passo primo PRIME_STEP=37 (coprimo con 1024) per selezionare le
 * coppie di campioni da ruotare:
 *
 *   coppia t: (i, (i + PRIME_STEP) % 1024)
 *
 * Poiché gcd(37, 1024) = 1, la sequenza di indici visita tutte le 1024
 * posizioni prima di ripetersi (teorema dell'orbita). Questo garantisce
 * una mescolanza GLOBALE dell'energia spettrale fin dal primo passo,
 * anziché una diffusione locale come nel Plugin 1.
 *
 * Effetto sonoro: cambiamento spettrale più ampio e immediato rispetto alla
 * catena locale. Utile per trasformazioni timbriche radicali o per esplorare
 * regioni distanti dello spazio SO(1024).
 *
 * SINTASSI CSOUND
 * ---------------
 *   aOut wtgivnlocal kfreq, kamp, itab,
 *                      kstart, kright, kleft, kboth,
 *                      ktheta0, kdecay, kmove
 *
 * PARAMETRI
 * ---------
 *   kfreq   (k) : frequenza [Hz]
 *   kamp    (k) : ampiezza
 *   itab    (i) : numero ftable (lunghezza >= 1024)
 *   kstart  (k) : indice di partenza [0..1023]
 *   kright  (k) : densità mixing locale forward (i,i+1) [0..1]
 *   kleft   (k) : densità mixing locale backward (i,i-1) [0..1]
 *   kboth   (k) : densità mixing non-locale (i, i+37) [0..1]
 *   ktheta0 (k) : angolo iniziale [radianti]
 *   kdecay  (k) : decadimento geometrico degli angoli [0..1]
 *   kmove   (k) : scaler globale trasformazione [0..1]
 *                 kmove=0 => wavetable base inalterata
 *                 kmove=1 => trasformazione completa
 *
 * NOTE
 * ----
 * - right/left/both controllano il NUMERO di rotazioni per k-block
 *   (0.0 => 0 rotazioni, 1.0 => MAX_STEPS=512 rotazioni).
 * - L'ordine applicazione è: both (non-locale) -> right -> left.
 * - MAX_STEPS=512 bilancia qualità sonora e budget CPU.
 * ============================================================================
 */

#include "csdl.h"
#include <math.h>
#include <string.h>

#define WT_LEN  1024
#define WT_MASK (WT_LEN - 1)

/* Passo primo coprimo con 1024: garantisce che le coppie non-locali
 * visitino tutte le 1024 posizioni prima di ripetersi. */
static const int PRIME_STEP = 37;

/* Budget massimo rotazioni per k-block per ciascuna direzione.
 * 512 è un buon compromesso tra ricchezza spettrale e CPU usage. */
static const int MAX_STEPS = 512;

/* ============================================================
 * Struttura dati dell'opcode
 * ============================================================ */
typedef struct {
  OPDS  h;

  MYFLT *out;

  /* parametri oscillatore */
  MYFLT *kfreq;
  MYFLT *kamp;
  MYFLT *itab;

  /* parametri dispersione */
  MYFLT *kstart;        /* indice di partenza [0..1023] */
  MYFLT *kright;        /* densità mixing locale forward [0..1] */
  MYFLT *kleft;         /* densità mixing locale backward [0..1] */
  MYFLT *kboth;         /* densità mixing non-locale [0..1] */
  MYFLT *ktheta0;       /* angolo iniziale [radianti] */
  MYFLT *kdecay;        /* decadimento geometrico [0..1] */
  MYFLT *kmove;         /* scaler globale [0..1] */

  /* stato interno */
  double phase;
  double sr;

  MYFLT wt[WT_LEN];
} WTGIVNLOCAL;


static inline int clampi(int v, int lo, int hi) {
  return (v < lo) ? lo : (v > hi) ? hi : v;
}

/* ============================================================
 * givens_inplace
 * x'[i] =  cos(t)*x[i] - sin(t)*x[j]
 * x'[j] =  sin(t)*x[i] + cos(t)*x[j]
 * ============================================================ */
static inline void givens_inplace(MYFLT *x, int i, int j, double theta)
{
  double c  = cos(theta), s = sin(theta);
  double xi = (double)x[i], xj = (double)x[j];
  x[i] = (MYFLT)(c*xi - s*xj);
  x[j] = (MYFLT)(s*xi + c*xj);
}

static void peak_normalize(MYFLT *x)
{
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
 * apply_right
 * steps rotazioni locali forward: (i, i+1), avanza i = j.
 * Effetto: diffusione spettrale locale in avanti.
 * ============================================================ */
static void apply_right(MYFLT *x, int start, int steps, double theta0, double decay)
{
  double th = theta0;
  int i = start & WT_MASK;
  for (int t = 0; t < steps; t++) {
    int j = (i + 1) & WT_MASK;
    givens_inplace(x, i, j, th);
    th *= decay;
    i = j;
  }
}

/* ============================================================
 * apply_left
 * steps rotazioni locali backward: (i, i-1), avanza i = j.
 * ============================================================ */
static void apply_left(MYFLT *x, int start, int steps, double theta0, double decay)
{
  double th = theta0;
  int i = start & WT_MASK;
  for (int t = 0; t < steps; t++) {
    int j = (i - 1) & WT_MASK;
    givens_inplace(x, i, j, th);
    th *= decay;
    i = j;
  }
}

/* ============================================================
 * apply_both_nonlocal
 * steps rotazioni NON-LOCALI a passo primo:
 *   coppia: (i, (i+PRIME_STEP) % 1024)
 *   drift: i avanza di 1 ad ogni passo (evita di ripetere lo stesso pattern)
 *
 * La combinazione passo-primo + drift garantisce che le rotazioni
 * raggiungano rapidamente tutto lo spazio della wavetable.
 * ============================================================ */
static void apply_both_nonlocal(MYFLT *x, int start, int steps, double theta0, double decay)
{
  double th = theta0;
  int i = start & WT_MASK;
  for (int t = 0; t < steps; t++) {
    int j = (i + PRIME_STEP) & WT_MASK;
    givens_inplace(x, i, j, th);
    th *= decay;
    i = (i + 1) & WT_MASK;  /* drift: avanza di 1 ad ogni passo */
  }
}

/* ============================================================
 * wtgivnlocal_init  [i-rate]
 * ============================================================ */
static int wtgivnlocal_init(CSOUND *csound, WTGIVNLOCAL *p)
{
  p->phase = 0.0;
  p->sr    = (double)csound->GetSr(csound);
  for (int i = 0; i < WT_LEN; i++) p->wt[i] = FL(0.0);
  return OK;
}

/* ============================================================
 * wtgivnlocal_perf  [a+k-rate]
 * Ogni k-block:
 *   1) Carica wavetable base dalla ftable
 *   2) Converte right/left/both in numero di rotazioni (0..MAX_STEPS)
 *   3) Applica: both_nonlocal -> right -> left
 *   4) Peak-normalizza
 *   5) Loop audio con interpolazione lineare
 * ============================================================ */
static int wtgivnlocal_perf(CSOUND *csound, WTGIVNLOCAL *p)
{
  MYFLT   *out   = p->out;
  uint32_t offset = p->h.insdshead->ksmps_offset;
  uint32_t early  = p->h.insdshead->ksmps_no_end;
  uint32_t nsmps  = CS_KSMPS;

  if (UNLIKELY(offset)) memset(out, 0, offset * sizeof(MYFLT));
  if (UNLIKELY(early)) { nsmps -= early; memset(&out[nsmps], 0, early * sizeof(MYFLT)); }

  /* 1) carica wavetable base */
  FUNC *ft = csound->FTnp2Find(csound, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) {
    for (uint32_t n = offset; n < nsmps; n++) out[n] = FL(0.0);
    return OK;
  }
  for (int i = 0; i < WT_LEN; i++) p->wt[i] = ft->ftable[i];

  /* 2) leggi e clamp parametri */
  int start = clampi((int)floor((double)*p->kstart + 0.5), 0, WT_LEN-1);

  double right  = (double)*p->kright;  if (right < 0) right = 0; if (right > 1) right = 1;
  double left   = (double)*p->kleft;   if (left  < 0) left  = 0; if (left  > 1) left  = 1;
  double both   = (double)*p->kboth;   if (both  < 0) both  = 0; if (both  > 1) both  = 1;
  double theta0 = (double)*p->ktheta0;
  double decay  = (double)*p->kdecay;  if (decay < 0) decay = 0; if (decay > 1) decay = 1;
  double move   = (double)*p->kmove;   if (move  < 0) move  = 0; if (move  > 1) move  = 1;

  /* kmove scala theta0: move=0 => nessuna rotazione */
  double theta_eff = theta0 * move;

  /* 3) converti densità in numero di rotazioni per k-block */
  int stepsR = (int)floor(right * (double)MAX_STEPS);
  int stepsL = (int)floor(left  * (double)MAX_STEPS);
  int stepsB = (int)floor(both  * (double)MAX_STEPS);

  /* 4) applica: non-locale prima (mescolanza globale), poi locali */
  if (stepsB > 0) apply_both_nonlocal(p->wt, start, stepsB, theta_eff, decay);
  if (stepsR > 0) apply_right        (p->wt, start, stepsR, theta_eff, decay);
  if (stepsL > 0) apply_left         (p->wt, start, stepsL, theta_eff, decay);

  /* 5) peak normalize */
  peak_normalize(p->wt);

  /* 6) loop audio con interpolazione lineare */
  double sr    = p->sr;
  double freq  = (double)*p->kfreq;  if (freq < 0) freq = 0;
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

    double y0 = (double)p->wt[i0];
    double y1 = (double)p->wt[i1];
    out[n] = (MYFLT)(amp * (y0 + (y1 - y0) * frac));
  }

  p->phase = phase;
  return OK;
}

/* ============================================================
 * Registrazione opcode
 * Firma: aOut wtgivnlocal kfreq, kamp, itab,
 *                           kstart, kright, kleft, kboth,
 *                           ktheta0, kdecay, kmove
 * ============================================================ */
static OENTRY localops[] = {
  {
    "wtgivnlocal",
    sizeof(WTGIVNLOCAL),
    0,
    3,
    "a",
    "kkikkkkkkk",
    (SUBR)wtgivnlocal_init,
    (SUBR)wtgivnlocal_perf
  }
};

LINKAGE
