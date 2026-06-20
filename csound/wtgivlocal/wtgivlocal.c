/*
 * wtgivlocal.c
 * ============================================================================
 * Csound plugin opcode: wtgivlocal
 *
 * DESCRIZIONE
 * -----------
 * Oscillatore wavetable con dispersione timbrica locale tramite ROTAZIONI DI
 * GIVENS su coppie adiacenti di campioni.
 *
 * Una rotazione di Givens G(i, j, theta) appartiene a SO(N): è una rotazione
 * nel piano (i,j) dello spazio R^N che preserva esattamente la norma L2
 * (energia) del vettore. Applicata a una wavetable di N=1024 campioni, mescola
 * solo i campioni in posizione i e j, lasciando invariati tutti gli altri:
 *
 *   x'[i] =  cos(theta)*x[i] - sin(theta)*x[j]
 *   x'[j] =  sin(theta)*x[i] + cos(theta)*x[j]
 *
 * Questo opcode applica una CATENA LOCALE: partendo da un indice (kstart),
 * la rotazione viene applicata a coppie adiacenti (k, k+1). L'angolo decade
 * geometricamente (theta_t = theta0 * decay^t * kmove), quindi i campioni
 * vicini a kstart sono più perturbati di quelli lontani.
 *
 * Effetto sonoro: diffusione spettrale lenta e controllata, modifica
 * preferenziale dei parziali legati alla regione kstart della wavetable.
 *
 * SINTASSI CSOUND
 * ---------------
 *   aOut wtgivlocal kfreq, kamp, itab,
 *                     kstart, kright, kleft, kboth,
 *                     ktheta0, kdecay, kmove
 *
 * PARAMETRI
 * ---------
 *   kfreq   (k) : frequenza di oscillazione [Hz]
 *   kamp    (k) : ampiezza di uscita
 *   itab    (i) : numero ftable Csound (lunghezza >= 1024 campioni)
 *   kstart  (k) : indice di partenza della catena [0..1022]
 *   kright  (k) : ampiezza catena verso destra [0..1]
 *   kleft   (k) : ampiezza catena verso sinistra [0..1]
 *   kboth   (k) : ampiezza catena bidirezionale [0..1]
 *   ktheta0 (k) : angolo iniziale [radianti]
 *   kdecay  (k) : decadimento geometrico degli angoli [0..1]
 *   kmove   (k) : scaler globale trasformazione [0..1]
 *                 kmove=0 => wavetable base inalterata
 *                 kmove=1 => trasformazione completa
 * ============================================================================
 */

#include "csdl.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define WT_LEN 1024

/* ============================================================
 * Struttura dati dell'opcode (una istanza per voce Csound)
 * ============================================================ */
typedef struct {
  OPDS  h;              /* header obbligatorio Csound */

  MYFLT *out;           /* buffer audio a-rate in uscita */

  /* parametri oscillatore */
  MYFLT *kfreq;         /* frequenza [Hz] */
  MYFLT *kamp;          /* ampiezza */
  MYFLT *itab;          /* numero ftable Csound */

  /* parametri dispersione */
  MYFLT *kstart;        /* indice di partenza catena [0..1022] */
  MYFLT *kright;        /* ampiezza catena destra [0..1] */
  MYFLT *kleft;         /* ampiezza catena sinistra [0..1] */
  MYFLT *kboth;         /* ampiezza catena bidirezionale [0..1] */
  MYFLT *ktheta0;       /* angolo iniziale [radianti] */
  MYFLT *kdecay;        /* decadimento geometrico [0..1] */
  MYFLT *kmove;         /* scaler globale [0..1]; 0=base, 1=pieno */

  /* stato interno */
  double phase;         /* fase oscillatore [0,1) */
  double sr;            /* sample rate [Hz] */

  MYFLT wt[WT_LEN];     /* copia di lavoro della wavetable trasformata */
} WTGIVLOCAL;


/* ============================================================
 * givens_inplace
 * Applica G(i,j,theta) al vettore x in-place.
 *
 *   x'[i] =  cos(theta)*x[i] - sin(theta)*x[j]
 *   x'[j] =  sin(theta)*x[i] + cos(theta)*x[j]
 *
 * Proprietà: conserva x[i]^2 + x[j]^2 (norma del piano).
 * ============================================================ */
static inline void givens_inplace(MYFLT *x, int i, int j, double theta)
{
  double c  = cos(theta);
  double s  = sin(theta);
  double xi = (double)x[i];
  double xj = (double)x[j];
  x[i] = (MYFLT)(c*xi - s*xj);
  x[j] = (MYFLT)(s*xi + c*xj);
}

static inline int clampi(int v, int lo, int hi)
{
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}

/* ============================================================
 * apply_chain_right
 * Catena locale verso destra a partire da start:
 *   (start,start+1), (start+1,start+2), ... , (WT_LEN-2, WT_LEN-1)
 * Angolo: theta_t = theta0 * decay^t
 * Effetto: diffusione spettrale che parte da kstart e si propaga a destra.
 * ============================================================ */
static void apply_chain_right(MYFLT *x, int start, double theta0, double decay)
{
  double th = theta0;
  for (int k = start; k <= WT_LEN - 2; k++) {
    givens_inplace(x, k, k+1, th);
    th *= decay;
  }
}

/* ============================================================
 * apply_chain_left
 * Prima applica una rotazione pivot (start, start+1), poi risale
 * verso sinistra: (start-1,start), (start-2,start-1), ..., (0,1).
 * ============================================================ */
static void apply_chain_left(MYFLT *x, int start, double theta0, double decay)
{
  double th = theta0;
  givens_inplace(x, start, start+1, th);
  th *= decay;
  for (int k = start; k >= 1; k--) {
    givens_inplace(x, k-1, k, th);
    th *= decay;
  }
}

/* ============================================================
 * apply_chain_both
 * Espansione bidirezionale simmetrica da start:
 *   t=0: (start,   start+1)
 *   t=1: (start-1, start)
 *   t=2: (start+1, start+2)
 *   t=3: (start-2, start-1)
 *   ...
 * Effetto: perturbazione centrata su kstart, decrescente verso i bordi.
 * ============================================================ */
static void apply_chain_both(MYFLT *x, int start, double theta0, double decay)
{
  double th = theta0;
  givens_inplace(x, start, start+1, th);
  th *= decay;

  int L = start - 1;
  int R = start + 1;
  while (L >= 0 || R <= WT_LEN - 2) {
    if (L >= 0) {
      givens_inplace(x, L, L+1, th);
      th *= decay;
      L -= 1;
    }
    if (R <= WT_LEN - 2) {
      givens_inplace(x, R, R+1, th);
      th *= decay;
      R += 1;
    }
  }
}

/* ============================================================
 * peak_normalize
 * Normalizza il vettore al picco assoluto = 1.0.
 * Necessaria perché le rotazioni conservano la norma L2 ma non il picco.
 * ============================================================ */
static void peak_normalize(MYFLT *x)
{
  double m = 0.0;
  for (int i = 0; i < WT_LEN; i++) {
    double a = fabs((double)x[i]);
    if (a > m) m = a;
  }
  if (m < 1e-12) return;
  double inv = 1.0 / m;
  for (int i = 0; i < WT_LEN; i++)
    x[i] = (MYFLT)((double)x[i] * inv);
}

/* ============================================================
 * wtgivlocal_init  [i-rate]
 * Inizializza fase e sample rate; azera la wavetable di lavoro.
 * ============================================================ */
static int wtgivlocal_init(CSOUND *csound, WTGIVLOCAL *p)
{
  p->phase = 0.0;
  p->sr    = (double)csound->GetSr(csound);
  for (int i = 0; i < WT_LEN; i++) p->wt[i] = FL(0.0);
  return OK;
}

/* ============================================================
 * wtgivlocal_perf  [a+k-rate]
 * Chiamata ogni k-block. Esegue:
 *   1) Carica wavetable base dalla ftable
 *   2) Applica catene Givens: both -> right -> left
 *   3) Peak-normalizza
 *   4) Loop audio con interpolazione lineare
 * ============================================================ */
static int wtgivlocal_perf(CSOUND *csound, WTGIVLOCAL *p)
{
  MYFLT   *out   = p->out;
  uint32_t offset = p->h.insdshead->ksmps_offset;
  uint32_t early  = p->h.insdshead->ksmps_no_end;
  uint32_t nsmps  = CS_KSMPS;

  if (UNLIKELY(offset))
    memset(out, 0, offset * sizeof(MYFLT));
  if (UNLIKELY(early)) {
    nsmps -= early;
    memset(&out[nsmps], 0, early * sizeof(MYFLT));
  }

  /* 1) carica wavetable base */
  FUNC *ft = csound->FTnp2Find(csound, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) {
    for (uint32_t n = offset; n < nsmps; n++) out[n] = FL(0.0);
    return OK;
  }
  for (int i = 0; i < WT_LEN; i++) p->wt[i] = ft->ftable[i];

  /* 2) leggi e clamp parametri */
  int start = clampi((int)floor((double)*p->kstart + 0.5), 0, WT_LEN - 2);

  double right  = (double)*p->kright;
  double left   = (double)*p->kleft;
  double both   = (double)*p->kboth;
  double theta0 = (double)*p->ktheta0;
  double decay  = (double)*p->kdecay;
  double move   = (double)*p->kmove;

  if (decay < 0.0) decay = 0.0;  if (decay > 1.0) decay = 1.0;
  if (move  < 0.0) move  = 0.0;  if (move  > 1.0) move  = 1.0;

  /* kmove scala theta0: move=0 => nessuna rotazione => wavetable base */
  double theta_eff = theta0 * move;

  /* 3) applica catene: both -> right -> left */
  if (both  != 0.0) apply_chain_both (p->wt, start, theta_eff * both,  decay);
  if (right != 0.0) apply_chain_right(p->wt, start, theta_eff * right, decay);
  if (left  != 0.0) apply_chain_left (p->wt, start, theta_eff * left,  decay);

  /* 4) peak normalize */
  peak_normalize(p->wt);

  /* 5) loop audio con interpolazione lineare */
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
    int    i1   = i0 + 1;
    if (i1 >= WT_LEN) i1 = 0;

    double y = (double)p->wt[i0] + ((double)p->wt[i1] - (double)p->wt[i0]) * frac;
    out[n] = (MYFLT)(amp * y);
  }

  p->phase = phase;
  return OK;
}

/* ============================================================
 * Registrazione opcode
 * Firma: aOut wtgivlocal kfreq, kamp, itab,
 *                          kstart, kright, kleft, kboth,
 *                          ktheta0, kdecay, kmove
 * ============================================================ */
static OENTRY localops[] = {
  {
    "wtgivlocal",
    sizeof(WTGIVLOCAL),
    0,
    3,
    "a",
    "kkikkkkkkk",
    (SUBR)wtgivlocal_init,
    (SUBR)wtgivlocal_perf
  }
};

LINKAGE
