/*
 * wtgivkura.c
 * ============================================================================
 * Csound plugin opcode: wtgivkura
 *
 * DESCRIPTION
 * -----------
 * Wavetable oscillator with continuous timbral navigation driven by the
 * KURAMOTO MODEL used as an angle source for Givens rotations.
 *
 * ── Kuramoto model ─────────────────────────────────────────────────────────
 * The Kuramoto model describes N coupled oscillators with phases theta_i:
 *
 *   dtheta_i/dt = omega_i + (K/N) * Σ_j sin(theta_j - theta_i)
 *
 * - omega_i : natural frequency of oscillator i (deterministic)
 * - K       : coupling strength (kboth)
 * - The sin(theta_j - theta_i) term tends to synchronize the phases
 *
 * For small K: each oscillator follows its own omega_i → controlled chaos
 * For large K: all oscillators synchronize → timbral convergence
 * The transition K_c = 2*sigma_omega / pi is the synchronization threshold
 *
 * ── Link to the wavetable ──────────────────────────────────────────────────
 * The theta[KDIM] vector (KDIM=64 oscillators) is used directly as the
 * angles for KDIM non-local Givens rotations on the wavetable.
 * Every k-block the Kuramoto system advances by dt = 0.002 * kmove,
 * producing a continuous trajectory in the space SO(1024).
 *
 * Key difference from Sobol:
 *   - Sobol is SPATIAL: same parameters → same timbre (static map)
 *   - Kuramoto is TEMPORAL: same parameters → same TRAJECTORY (flow)
 *
 * CSOUND SYNTAX
 * -------------
 *   aOut wtgivkura kfreq, kamp, itab,
 *                       kstart, kright, kleft, kboth,
 *                       ktheta0, kdecay, kmove
 *
 * PARAMETERS
 * ----------
 *   kfreq   (k) : frequency [Hz]
 *   kamp    (k) : amplitude
 *   itab    (i) : ftable number (length >= 1024)
 *   kstart  (k) : starting index in the wavetable [0..1023]
 *   kright  (k) : rightward directional bias of the omega profile [0..1]
 *   kleft   (k) : leftward directional bias of the omega profile [0..1]
 *   kboth   (k) : Kuramoto coupling strength K [0..10]
 *   ktheta0 (k) : scale of the natural frequencies omega [radians]
 *   kdecay  (k) : decay of the omega profile along the indices [0..1]
 *   kmove   (k) : speed of temporal evolution [0..10]
 *                 kmove=0 => Kuramoto system stationary => static timbre
 *
 * NOTES
 * -----
 * - The theta[] phases persist across k-blocks: the timbre evolves over time.
 * - Integration uses RK2 (second-order Runge-Kutta) for stability.
 * - KDIM=64 balances richness of the flow against CPU budget.
 * ============================================================================
 */

#include "csdl.h"
#include <math.h>
#include <string.h>

#define WT_LEN     1024
#define WT_MASK    (WT_LEN - 1)
#define KDIM       64      /* number of Kuramoto oscillators and Givens rotations */
#define PRIME_STEP 37      /* non-local step, coprime with 1024 */

static inline double clampd(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}
static inline int clampi(int v, int lo, int hi) {
  return v < lo ? lo : (v > hi ? hi : v);
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
 * kuramoto_rhs
 * Computes the derivatives of the Kuramoto system:
 *   dtheta_i = omega_i + (K/N) * Σ_j sin(theta_j - theta_i)
 *
 * The coupling term is the average of the sines of the phase
 * differences, which tends to align the phases when K is large.
 * ============================================================ */
static void kuramoto_rhs(const double *theta, const double *omega,
                          double K, double *dtheta)
{
  for (int i = 0; i < KDIM; i++) {
    double sum = 0.0;
    double ti  = theta[i];
    for (int j = 0; j < KDIM; j++)
      sum += sin(theta[j] - ti);
    dtheta[i] = omega[i] + K * (sum / (double)KDIM);
  }
}

/* ============================================================
 * kuramoto_step_rk2
 * RK2 integration (midpoint method) for numerical stability:
 *   k1 = f(theta)
 *   k2 = f(theta + 0.5*dt*k1)   ← evaluated at the midpoint
 *   theta += dt * k2
 *
 * RK2 is more stable than Euler for oscillatory systems.
 * ============================================================ */
static void kuramoto_step_rk2(double *theta, const double *omega,
                               double K, double dt)
{
  double k1[KDIM], k2[KDIM], mid[KDIM];
  kuramoto_rhs(theta, omega, K, k1);
  for (int i = 0; i < KDIM; i++) mid[i] = theta[i] + 0.5*dt*k1[i];
  kuramoto_rhs(mid, omega, K, k2);
  for (int i = 0; i < KDIM; i++) theta[i] += dt * k2[i];
}

/* ============================================================
 * make_omega
 * Builds the natural frequency profile omega[KDIM].
 *
 * Base profile: omega_i = theta0 * decay^i
 *   → geometric decay along the indices (like the local chains)
 *
 * Directional bias tilt = 0.25*(right-left):
 *   omega_i *= (1 + tilt * pos_i)   where pos_i ∈ [-1, +1]
 *   → right>left: higher-index oscillators get a larger omega
 *   → left>right: lower-index oscillators get a larger omega
 *
 * right and left influence the "direction" of the Kuramoto flow,
 * not the number of rotations as in the other opcodes.
 * ============================================================ */
static void make_omega(double *omega, double theta0, double decay,
                       double right, double left)
{
  double tilt = 0.25 * (right - left);
  for (int i = 0; i < KDIM; i++) {
    double base = theta0 * pow(decay, (double)i);
    double pos  = -1.0 + 2.0 * (double)i / (double)(KDIM - 1);
    omega[i] = base * (1.0 + tilt * pos);
  }
}

/* ============================================================
 * Opcode data structure
 * ============================================================ */
typedef struct {
  OPDS h;
  MYFLT *out;

  MYFLT *kfreq, *kamp, *itab;
  MYFLT *kstart, *kright, *kleft, *kboth, *ktheta0, *kdecay, *kmove;

  double phase, sr;

  MYFLT base[WT_LEN];     /* original wavetable */
  MYFLT wt[WT_LEN];       /* transformed wavetable */

  double theta[KDIM];     /* Kuramoto phases — PERSIST across k-blocks!
                           * This is the "memory" of the timbral flow. */
} OP;

static int load_base(CSOUND *cs, OP *p) {
  FUNC *ft = cs->FTnp2Find(cs, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) return NOTOK;
  for (int i = 0; i < WT_LEN; i++) p->base[i] = ft->ftable[i];
  return OK;
}

/* ============================================================
 * build_wavetable
 * Uses the current KDIM Kuramoto phases as angles for KDIM
 * non-local Givens rotations (prime step) on the wavetable.
 *
 * i starts at start, advancing by 1 (drift) on each rotation.
 * j = (i + PRIME_STEP) % 1024
 * ============================================================ */
static void build_wavetable(OP *p, int start) {
  memcpy(p->wt, p->base, WT_LEN * sizeof(MYFLT));
  int i = start & WT_MASK;
  for (int k = 0; k < KDIM; k++) {
    int j = (i + PRIME_STEP) & WT_MASK;
    givens_inplace(p->wt, i, j, p->theta[k]);
    i = (i + 1) & WT_MASK;
  }
  peak_normalize(p->wt);
}

/* Audio loop with linear interpolation */
static void render(OP *p, MYFLT *out, uint32_t offset, uint32_t nsmps) {
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

    out[n] = (MYFLT)(amp * ((double)p->wt[i0] +
             ((double)p->wt[i1] - (double)p->wt[i0]) * frac));
  }
  p->phase = phase;
}

/* ============================================================
 * op_init  [i-rate]
 * Initializes Kuramoto phases to zero: theta_i=0 => base wavetable.
 * ============================================================ */
static int op_init(CSOUND *cs, OP *p) {
  p->phase = 0.0;
  p->sr    = (double)cs->GetSr(cs);
  for (int i = 0; i < KDIM; i++) p->theta[i] = 0.0;
  return OK;
}

/* ============================================================
 * op_perf  [a+k-rate]
 * Every k-block:
 *   1) Load the base wavetable
 *   2) Build the omega profile from the parameters
 *   3) Advance the Kuramoto system by dt = 0.002 * kmove
 *   4) Build the wavetable with the current phases
 *   5) Render audio
 * ============================================================ */
static int op_perf(CSOUND *cs, OP *p) {
  MYFLT   *out   = p->out;
  uint32_t offset = p->h.insdshead->ksmps_offset;
  uint32_t early  = p->h.insdshead->ksmps_no_end;
  uint32_t nsmps  = CS_KSMPS;

  if (UNLIKELY(offset)) memset(out, 0, offset * sizeof(MYFLT));
  if (UNLIKELY(early))  { nsmps -= early; memset(&out[nsmps], 0, early * sizeof(MYFLT)); }

  if (UNLIKELY(load_base(cs, p) != OK)) {
    for (uint32_t n = offset; n < nsmps; n++) out[n] = FL(0.0);
    return OK;
  }

  int start = clampi((int)floor((double)*p->kstart + 0.5), 0, WT_LEN-1);

  double right  = clampd((double)*p->kright,  0.0, 1.0);
  double left   = clampd((double)*p->kleft,   0.0, 1.0);
  double K      = clampd((double)*p->kboth,   0.0, 10.0);  /* coupling strength */
  double theta0 = (double)*p->ktheta0;                     /* omega scale */
  double decay  = clampd((double)*p->kdecay,  0.0, 1.0);
  double move   = clampd((double)*p->kmove,   0.0, 10.0);  /* dt speed */

  /* build the omega profile */
  double omega[KDIM];
  make_omega(omega, theta0, decay, right, left);

  /* advance the Kuramoto system:
   * small dt = slow, stable evolution
   * kmove=0 => dt=0 => theta stationary => timbre unchanged */
  double dt = 0.002 * move;
  kuramoto_step_rk2(p->theta, omega, K, dt);

  /* build the wavetable with the updated phases */
  build_wavetable(p, start);
  render(p, out, offset, nsmps);
  return OK;
}

/* ============================================================
 * Opcode registration
 * Signature: aOut wtgivkura kfreq, kamp, itab,
 *                            kstart, kright, kleft, kboth,
 *                            ktheta0, kdecay, kmove
 * ============================================================ */
static OENTRY localops[] = {
  { "wtgivkura", sizeof(OP), 0, 3, "a", "kkikkkkkkk",
    (SUBR)op_init, (SUBR)op_perf }
};

LINKAGE