/*
 * wtgivlocal.c
 * ============================================================================
 * Csound plugin opcode: wtgivlocal
 *
 * DESCRIPTION
 * -----------
 * Wavetable oscillator with local timbral dispersion via GIVENS ROTATIONS
 * on adjacent pairs of samples.
 *
 * A Givens rotation G(i, j, theta) belongs to SO(N): it is a rotation
 * in the (i,j) plane of R^N space that exactly preserves the L2 norm
 * (energy) of the vector. Applied to a wavetable of N=1024 samples, it
 * mixes only the samples at positions i and j, leaving all others unchanged:
 *
 *   x'[i] =  cos(theta)*x[i] - sin(theta)*x[j]
 *   x'[j] =  sin(theta)*x[i] + cos(theta)*x[j]
 *
 * This opcode applies a LOCAL CHAIN: starting from an index (kstart),
 * the rotation is applied to adjacent pairs (k, k+1). The angle decays
 * geometrically (theta_t = theta0 * decay^t * kmove), so samples near
 * kstart are perturbed more than those farther away.
 *
 * Sonic effect: slow, controlled spectral diffusion, with preferential
 * modification of the partials tied to the kstart region of the wavetable.
 *
 * CSOUND SYNTAX
 * -------------
 *   aOut wtgivlocal kfreq, kamp, itab,
 *                     kstart, kright, kleft, kboth,
 *                     ktheta0, kdecay, kmove
 *
 * PARAMETERS
 * ----------
 *   kfreq   (k) : oscillation frequency [Hz]
 *   kamp    (k) : output amplitude
 *   itab    (i) : Csound ftable number (length >= 1024 samples)
 *   kstart  (k) : chain starting index [0..1022]
 *   kright  (k) : rightward chain amplitude [0..1]
 *   kleft   (k) : leftward chain amplitude [0..1]
 *   kboth   (k) : bidirectional chain amplitude [0..1]
 *   ktheta0 (k) : initial angle [radians]
 *   kdecay  (k) : geometric decay of the angles [0..1]
 *   kmove   (k) : global transformation scaler [0..1]
 *                 kmove=0 => base wavetable unaltered
 *                 kmove=1 => full transformation
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
 * Opcode data structure (one instance per Csound voice)
 * ============================================================ */
typedef struct {
  OPDS  h;              /* mandatory Csound header */

  MYFLT *out;           /* a-rate audio output buffer */

  /* oscillator parameters */
  MYFLT *kfreq;         /* frequency [Hz] */
  MYFLT *kamp;          /* amplitude */
  MYFLT *itab;          /* Csound ftable number */

  /* dispersion parameters */
  MYFLT *kstart;        /* chain starting index [0..1022] */
  MYFLT *kright;        /* right chain amplitude [0..1] */
  MYFLT *kleft;         /* left chain amplitude [0..1] */
  MYFLT *kboth;         /* bidirectional chain amplitude [0..1] */
  MYFLT *ktheta0;       /* initial angle [radians] */
  MYFLT *kdecay;        /* geometric decay [0..1] */
  MYFLT *kmove;         /* global scaler [0..1]; 0=base, 1=full */

  /* internal state */
  double phase;         /* oscillator phase [0,1) */
  double sr;            /* sample rate [Hz] */

  MYFLT wt[WT_LEN];     /* working copy of the transformed wavetable */
} WTGIVLOCAL;


/* ============================================================
 * givens_inplace
 * Applies G(i,j,theta) to the vector x in-place.
 *
 *   x'[i] =  cos(theta)*x[i] - sin(theta)*x[j]
 *   x'[j] =  sin(theta)*x[i] + cos(theta)*x[j]
 *
 * Property: preserves x[i]^2 + x[j]^2 (plane norm).
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
 * Local chain moving rightward starting from start:
 *   (start,start+1), (start+1,start+2), ... , (WT_LEN-2, WT_LEN-1)
 * Angle: theta_t = theta0 * decay^t
 * Effect: spectral diffusion starting at kstart and propagating rightward.
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
 * First applies a pivot rotation (start, start+1), then moves back
 * leftward: (start-1,start), (start-2,start-1), ..., (0,1).
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
 * Symmetric bidirectional expansion from start:
 *   t=0: (start,   start+1)
 *   t=1: (start-1, start)
 *   t=2: (start+1, start+2)
 *   t=3: (start-2, start-1)
 *   ...
 * Effect: perturbation centered on kstart, decreasing toward the edges.
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
 * Normalizes the vector to an absolute peak of 1.0.
 * Necessary because the rotations preserve the L2 norm but not the peak.
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
 * Initializes phase and sample rate; zeroes the working wavetable.
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
 * Called every k-block. Performs:
 *   1) Load the base wavetable from the ftable
 *   2) Apply Givens chains: both -> right -> left
 *   3) Peak-normalize
 *   4) Audio loop with linear interpolation
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

  /* 1) load base wavetable */
  FUNC *ft = csound->FTnp2Find(csound, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) {
    for (uint32_t n = offset; n < nsmps; n++) out[n] = FL(0.0);
    return OK;
  }
  for (int i = 0; i < WT_LEN; i++) p->wt[i] = ft->ftable[i];

  /* 2) read and clamp parameters */
  int start = clampi((int)floor((double)*p->kstart + 0.5), 0, WT_LEN - 2);

  double right  = (double)*p->kright;
  double left   = (double)*p->kleft;
  double both   = (double)*p->kboth;
  double theta0 = (double)*p->ktheta0;
  double decay  = (double)*p->kdecay;
  double move   = (double)*p->kmove;

  if (decay < 0.0) decay = 0.0;  if (decay > 1.0) decay = 1.0;
  if (move  < 0.0) move  = 0.0;  if (move  > 1.0) move  = 1.0;

  /* kmove scales theta0: move=0 => no rotation => base wavetable */
  double theta_eff = theta0 * move;

  /* 3) apply chains: both -> right -> left */
  if (both  != 0.0) apply_chain_both (p->wt, start, theta_eff * both,  decay);
  if (right != 0.0) apply_chain_right(p->wt, start, theta_eff * right, decay);
  if (left  != 0.0) apply_chain_left (p->wt, start, theta_eff * left,  decay);

  /* 4) peak normalize */
  peak_normalize(p->wt);

  /* 5) audio loop with linear interpolation */
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
 * Opcode registration
 * Signature: aOut wtgivlocal kfreq, kamp, itab,
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