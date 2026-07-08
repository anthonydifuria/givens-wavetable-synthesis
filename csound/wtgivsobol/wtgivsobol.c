/*
 * wtgivsobol.c
 * ============================================================================
 * Csound plugin opcode: wtgivsobol
 *
 * DESCRIPTION
 * -----------
 * Wavetable oscillator with timbral mixing driven by SOBOL SEQUENCES
 * applied as angles for Givens rotations.
 *
 * ── What is a Sobol sequence? ──────────────────────────────────────────────
 * A Sobol sequence is a low-discrepancy quasi-random sequence: the points
 * are "more uniformly distributed" than purely pseudo-random points.
 * In 1D, the Sobol sequence fills [0,1) very evenly, avoiding clusters
 * and gaps. This produces rotation angles distributed in a nearly
 * uniform way over [-1,1], ensuring that every region of the SO(1024)
 * timbral space gets explored in a balanced manner.
 *
 * 1D Sobol algorithm (Gray-code increment):
 *   x_n = x_{n-1} XOR v_c     where c = trailing_zeros(n)
 *   v_c = 1 << (31 - c)        (simplified direction numbers)
 *
 * ── DETERMINISTIC behavior ─────────────────────────────────────────────────
 * IMPORTANT: the Sobol generator is RESET at every k-block (seed=0).
 * This means the same parameters ALWAYS produce the same wavetable —
 * identical behavior to the Python version.
 * It is a stable map from parameters → timbre, not a random drift.
 * kmove=0 returns exactly the base wavetable.
 *
 * CSOUND SYNTAX
 * -------------
 *   aOut wtgivsobol kfreq, kamp, itab,
 *                       kstart, kright, kleft, kboth,
 *                       ktheta0, kdecay, kmove
 *
 * PARAMETERS
 * ----------
 *   kfreq   (k) : frequency [Hz]
 *   kamp    (k) : amplitude
 *   itab    (i) : ftable number (length >= 1024)
 *   kstart  (k) : starting index [0..1023]
 *   kright  (k) : forward local mixing density [0..1]
 *   kleft   (k) : backward local mixing density [0..1]
 *   kboth   (k) : non-local mixing density (prime step) [0..1]
 *   ktheta0 (k) : angle scale [radians]
 *   kdecay  (k) : geometric decay [0..1]
 *   kmove   (k) : transformation intensity [0..1]
 *                 kmove=0 => base wavetable unaltered
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
 * SOBOL 1D — deterministic low-discrepancy generator
 *
 * State: integer counter i + accumulator x (32-bit XOR).
 * On each call to sobol_next01():
 *   1) increment i
 *   2) compute c = number of trailing zeros of i (count trailing zeros)
 *   3) direction number: v = 1 << (31 - c)
 *   4) x ^= v
 *   5) return x / 2^32  in [0,1)
 *
 * The resulting sequence has discrepancy O(log(N)/N), much better
 * than pseudo-random O(1/sqrt(N)).
 * ============================================================ */
typedef struct {
  uint32_t i;   /* call counter (starts at 0) */
  uint32_t x;   /* XOR accumulator (current state) */
} SOBOL1D;

/* Counts trailing zero bits of n>0 */
static inline int ctz32(uint32_t n) {
  int c = 0;
  while ((n & 1u) == 0u) { n >>= 1; c++; }
  return c;
}

/* Next sample in [0,1) */
static inline double sobol_next01(SOBOL1D *s) {
  s->i += 1u;
  int c = ctz32(s->i);
  uint32_t v = (uint32_t)(1u << (31 - c));
  s->x ^= v;
  return (double)s->x / 4294967296.0;   /* / 2^32 */
}

/* Deterministic reset: same seed => same sequence */
static inline void sobol_reset(SOBOL1D *s) {
  s->i = 0u;
  s->x = 0u;
}

/* ============================================================
 * Opcode data structure
 * ============================================================ */
typedef struct {
  OPDS h;
  MYFLT *out;

  MYFLT *kfreq, *kamp, *itab;
  MYFLT *kstart, *kright, *kleft, *kboth, *ktheta0, *kdecay, *kmove;

  double phase;
  double sr;

  MYFLT base[WT_LEN];   /* static copy of the original wavetable */
  MYFLT wt[WT_LEN];     /* transformed working copy */

  SOBOL1D sob;           /* Sobol generator (reset every k-block) */
} OP;

/* Loads the base wavetable from the Csound ftable */
static int load_base(CSOUND *csound, OP *p) {
  FUNC *ft = csound->FTnp2Find(csound, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) return NOTOK;
  for (int i = 0; i < WT_LEN; i++) p->base[i] = ft->ftable[i];
  return OK;
}

/* ============================================================
 * next_theta
 * Generates the next rotation angle from the Sobol sequence.
 *
 * u  ∈ [0,1)  from Sobol
 * z  = 2u-1  ∈ [-1,1]  (centered on zero)
 * theta = move * theta0 * decay^t * z
 *
 * The geometric decay decay^t makes the first rotations (near
 * kstart) have larger angles than the later ones.
 * ============================================================ */
static inline double next_theta(OP *p, int t, double move, double theta0, double decay) {
  double u = sobol_next01(&p->sob);
  double z = 2.0*u - 1.0;
  return move * theta0 * pow(decay, (double)t) * z;
}

/* ============================================================
 * build_wavetable
 * Builds the transformed wavetable by applying Givens rotations
 * with angles driven by the Sobol sequence.
 *
 * DETERMINISTIC RESET: sobol_reset() is called at the start,
 * so the same combination of parameters always produces the
 * same wavetable (identical behavior to the Python code).
 * ============================================================ */
static void build_wavetable(OP *p, int start,
                            double right, double left, double both,
                            double theta0, double decay, double move)
{
  /* Sobol reset: ESSENTIAL for deterministic behavior.
   * Without the reset, the Sobol generator would keep advancing
   * from one k-block to the next, producing uncontrolled drift
   * even with fixed parameters. */
  sobol_reset(&p->sob);

  memcpy(p->wt, p->base, WT_LEN * sizeof(MYFLT));

  const int PRIME_STEP = 37;   /* coprime with 1024: guarantees global visitation */
  const int MAX_STEPS  = 512;  /* CPU budget per k-block */

  int stepsB = (int)floor(both  * (double)MAX_STEPS);
  int stepsR = (int)floor(right * (double)MAX_STEPS);
  int stepsL = (int)floor(left  * (double)MAX_STEPS);

  int t = 0;   /* global step counter (used for the decay) */

  /* BOTH: non-local prime-step mixing */
  int i = start & WT_MASK;
  for (int s = 0; s < stepsB; s++) {
    int j = (i + PRIME_STEP) & WT_MASK;
    double th = next_theta(p, t++, move, theta0, decay);
    givens_inplace(p->wt, i, j, th);
    i = (i + 1) & WT_MASK;   /* drift to avoid repeating the same pairs */
  }

  /* RIGHT: forward local mixing */
  i = start & WT_MASK;
  for (int s = 0; s < stepsR; s++) {
    int j = (i + 1) & WT_MASK;
    double th = next_theta(p, t++, move, theta0, decay);
    givens_inplace(p->wt, i, j, th);
    i = j;
  }

  /* LEFT: backward local mixing */
  i = start & WT_MASK;
  for (int s = 0; s < stepsL; s++) {
    int j = (i - 1) & WT_MASK;
    double th = next_theta(p, t++, move, theta0, decay);
    givens_inplace(p->wt, i, j, th);
    i = j;
  }

  peak_normalize(p->wt);
}

/* Audio loop with linear interpolation */
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
  sobol_reset(&p->sob);   /* deterministic seed = 0 */
  return OK;
}

/* ============================================================
 * op_perf  [a+k-rate]
 * Every k-block:
 *   1) Load the base wavetable
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
 * Opcode registration
 * Signature: aOut wtgivsobol kfreq, kamp, itab,
 *                            kstart, kright, kleft, kboth,
 *                            ktheta0, kdecay, kmove
 * ============================================================ */
static OENTRY localops[] = {
  { "wtgivsobol", sizeof(OP), 0, 3, "a", "kkikkkkkkk",
    (SUBR)op_init, (SUBR)op_perf }
};

LINKAGE