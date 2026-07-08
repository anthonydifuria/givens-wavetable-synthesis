/*
 * wtgivnlocal.c
 * ============================================================================
 * Csound plugin opcode: wtgivnlocal
 *
 * DESCRIPTION
 * -----------
 * Wavetable oscillator with NON-LOCAL timbral mixing via prime-step
 * GIVENS ROTATIONS.
 *
 * Unlike wtgivlocal (which operates on adjacent pairs), this opcode
 * uses a prime step PRIME_STEP=37 (coprime with 1024) to select the
 * pairs of samples to rotate:
 *
 *   pair t: (i, (i + PRIME_STEP) % 1024)
 *
 * Since gcd(37, 1024) = 1, the sequence of indices visits all 1024
 * positions before repeating (orbit theorem). This guarantees a GLOBAL
 * mixing of spectral energy from the very first step, rather than the
 * local diffusion seen in Plugin 1.
 *
 * Sonic effect: broader and more immediate spectral change compared to
 * the local chain. Useful for radical timbral transformations or for
 * exploring distant regions of the SO(1024) space.
 *
 * CSOUND SYNTAX
 * -------------
 *   aOut wtgivnlocal kfreq, kamp, itab,
 *                      kstart, kright, kleft, kboth,
 *                      ktheta0, kdecay, kmove
 *
 * PARAMETERS
 * ----------
 *   kfreq   (k) : frequency [Hz]
 *   kamp    (k) : amplitude
 *   itab    (i) : ftable number (length >= 1024)
 *   kstart  (k) : starting index [0..1023]
 *   kright  (k) : forward local mixing density (i,i+1) [0..1]
 *   kleft   (k) : backward local mixing density (i,i-1) [0..1]
 *   kboth   (k) : non-local mixing density (i, i+37) [0..1]
 *   ktheta0 (k) : initial angle [radians]
 *   kdecay  (k) : geometric decay of the angles [0..1]
 *   kmove   (k) : global transformation scaler [0..1]
 *                 kmove=0 => base wavetable unaltered
 *                 kmove=1 => full transformation
 *
 * NOTES
 * -----
 * - right/left/both control the NUMBER of rotations per k-block
 *   (0.0 => 0 rotations, 1.0 => MAX_STEPS=512 rotations).
 * - Application order is: both (non-local) -> right -> left.
 * - MAX_STEPS=512 balances sound quality and CPU budget.
 * ============================================================================
 */

#include "csdl.h"
#include <math.h>
#include <string.h>

#define WT_LEN  1024
#define WT_MASK (WT_LEN - 1)

/* Prime step coprime with 1024: guarantees that the non-local pairs
 * visit all 1024 positions before repeating. */
static const int PRIME_STEP = 37;

/* Maximum rotation budget per k-block for each direction.
 * 512 is a good compromise between spectral richness and CPU usage. */
static const int MAX_STEPS = 512;

/* ============================================================
 * Opcode data structure
 * ============================================================ */
typedef struct {
  OPDS  h;

  MYFLT *out;

  /* oscillator parameters */
  MYFLT *kfreq;
  MYFLT *kamp;
  MYFLT *itab;

  /* dispersion parameters */
  MYFLT *kstart;        /* starting index [0..1023] */
  MYFLT *kright;        /* forward local mixing density [0..1] */
  MYFLT *kleft;         /* backward local mixing density [0..1] */
  MYFLT *kboth;         /* non-local mixing density [0..1] */
  MYFLT *ktheta0;       /* initial angle [radians] */
  MYFLT *kdecay;        /* geometric decay [0..1] */
  MYFLT *kmove;         /* global scaler [0..1] */

  /* internal state */
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
 * steps forward local rotations: (i, i+1), advance i = j.
 * Effect: local spectral diffusion moving forward.
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
 * steps backward local rotations: (i, i-1), advance i = j.
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
 * steps NON-LOCAL prime-step rotations:
 *   pair: (i, (i+PRIME_STEP) % 1024)
 *   drift: i advances by 1 at each step (avoids repeating the same pattern)
 *
 * The combination of prime-step + drift guarantees that the rotations
 * quickly reach the entire span of the wavetable.
 * ============================================================ */
static void apply_both_nonlocal(MYFLT *x, int start, int steps, double theta0, double decay)
{
  double th = theta0;
  int i = start & WT_MASK;
  for (int t = 0; t < steps; t++) {
    int j = (i + PRIME_STEP) & WT_MASK;
    givens_inplace(x, i, j, th);
    th *= decay;
    i = (i + 1) & WT_MASK;  /* drift: advances by 1 at each step */
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
 * Every k-block:
 *   1) Load the base wavetable from the ftable
 *   2) Convert right/left/both into a number of rotations (0..MAX_STEPS)
 *   3) Apply: both_nonlocal -> right -> left
 *   4) Peak-normalize
 *   5) Audio loop with linear interpolation
 * ============================================================ */
static int wtgivnlocal_perf(CSOUND *csound, WTGIVNLOCAL *p)
{
  MYFLT   *out   = p->out;
  uint32_t offset = p->h.insdshead->ksmps_offset;
  uint32_t early  = p->h.insdshead->ksmps_no_end;
  uint32_t nsmps  = CS_KSMPS;

  if (UNLIKELY(offset)) memset(out, 0, offset * sizeof(MYFLT));
  if (UNLIKELY(early)) { nsmps -= early; memset(&out[nsmps], 0, early * sizeof(MYFLT)); }

  /* 1) load base wavetable */
  FUNC *ft = csound->FTnp2Find(csound, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) {
    for (uint32_t n = offset; n < nsmps; n++) out[n] = FL(0.0);
    return OK;
  }
  for (int i = 0; i < WT_LEN; i++) p->wt[i] = ft->ftable[i];

  /* 2) read and clamp parameters */
  int start = clampi((int)floor((double)*p->kstart + 0.5), 0, WT_LEN-1);

  double right  = (double)*p->kright;  if (right < 0) right = 0; if (right > 1) right = 1;
  double left   = (double)*p->kleft;   if (left  < 0) left  = 0; if (left  > 1) left  = 1;
  double both   = (double)*p->kboth;   if (both  < 0) both  = 0; if (both  > 1) both  = 1;
  double theta0 = (double)*p->ktheta0;
  double decay  = (double)*p->kdecay;  if (decay < 0) decay = 0; if (decay > 1) decay = 1;
  double move   = (double)*p->kmove;   if (move  < 0) move  = 0; if (move  > 1) move  = 1;

  /* kmove scales theta0: move=0 => no rotation */
  double theta_eff = theta0 * move;

  /* 3) convert density into number of rotations per k-block */
  int stepsR = (int)floor(right * (double)MAX_STEPS);
  int stepsL = (int)floor(left  * (double)MAX_STEPS);
  int stepsB = (int)floor(both  * (double)MAX_STEPS);

  /* 4) apply: non-local first (global mixing), then local */
  if (stepsB > 0) apply_both_nonlocal(p->wt, start, stepsB, theta_eff, decay);
  if (stepsR > 0) apply_right        (p->wt, start, stepsR, theta_eff, decay);
  if (stepsL > 0) apply_left         (p->wt, start, stepsL, theta_eff, decay);

  /* 5) peak normalize */
  peak_normalize(p->wt);

  /* 6) audio loop with linear interpolation */
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
 * Opcode registration
 * Signature: aOut wtgivnlocal kfreq, kamp, itab,
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