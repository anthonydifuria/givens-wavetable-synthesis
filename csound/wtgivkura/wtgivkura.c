/*
 * wtgivkura.c
 * ============================================================================
 * Csound plugin opcode: wtgivkura
 *
 * DESCRIZIONE
 * -----------
 * Oscillatore wavetable con navigazione timbrica continua guidata dal MODELLO
 * DI KURAMOTO applicato come sorgente di angoli per rotazioni di Givens.
 *
 * ── Modello di Kuramoto ───────────────────────────────────────────────────
 * Il modello di Kuramoto descrive N oscillatori accoppiati con fasi theta_i:
 *
 *   dtheta_i/dt = omega_i + (K/N) * Σ_j sin(theta_j - theta_i)
 *
 * - omega_i : frequenza naturale dell'oscillatore i (deterministica)
 * - K       : forza di accoppiamento (kboth)
 * - Il termine sin(theta_j - theta_i) tende a sincronizzare le fasi
 *
 * Per K piccolo: ogni oscillatore segue la sua omega_i → caos controllato
 * Per K grande: tutti si sincronizzano → convergenza timbrica
 * La transizione K_c = 2*sigma_omega / pi è la soglia di sincronizzazione
 *
 * ── Collegamento con la wavetable ────────────────────────────────────────
 * Il vettore theta[KDIM] (KDIM=64 oscillatori) viene usato direttamente
 * come angoli per KDIM rotazioni di Givens non-locali sulla wavetable.
 * Ogni k-block il sistema Kuramoto avanza di dt = 0.002 * kmove,
 * producendo una traiettoria continua nello spazio SO(1024).
 *
 * Differenza fondamentale rispetto a Sobol:
 *   - Sobol è SPAZIALE: stessi parametri → stesso timbro (mappa statica)
 *   - Kuramoto è TEMPORALE: stessi parametri → stessa TRAIETTORIA (flusso)
 *
 * SINTASSI CSOUND
 * ---------------
 *   aOut wtgivkura kfreq, kamp, itab,
 *                       kstart, kright, kleft, kboth,
 *                       ktheta0, kdecay, kmove
 *
 * PARAMETRI
 * ---------
 *   kfreq   (k) : frequenza [Hz]
 *   kamp    (k) : ampiezza
 *   itab    (i) : numero ftable (lunghezza >= 1024)
 *   kstart  (k) : indice di partenza nella wavetable [0..1023]
 *   kright  (k) : bias direzionale destra del profilo omega [0..1]
 *   kleft   (k) : bias direzionale sinistra del profilo omega [0..1]
 *   kboth   (k) : forza di accoppiamento Kuramoto K [0..10]
 *   ktheta0 (k) : scala delle frequenze naturali omega [radianti]
 *   kdecay  (k) : decadimento del profilo omega lungo gli indici [0..1]
 *   kmove   (k) : velocità di evoluzione temporale [0..10]
 *                 kmove=0 => sistema Kuramoto fermo => timbro statico
 *
 * NOTE
 * ----
 * - Le fasi theta[] persistono tra i k-block: il timbro evolve nel tempo.
 * - L'integrazione usa RK2 (Runge-Kutta secondo ordine) per stabilità.
 * - KDIM=64 bilancia ricchezza del flusso e budget CPU.
 * ============================================================================
 */

#include "csdl.h"
#include <math.h>
#include <string.h>

#define WT_LEN     1024
#define WT_MASK    (WT_LEN - 1)
#define KDIM       64      /* numero oscillatori Kuramoto e rotazioni Givens */
#define PRIME_STEP 37      /* passo non-locale, coprimo con 1024 */

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
 * Calcola le derivate del sistema Kuramoto:
 *   dtheta_i = omega_i + (K/N) * Σ_j sin(theta_j - theta_i)
 *
 * Il termine di accoppiamento è la media dei sin delle differenze
 * di fase, che tende ad allineare le fasi quando K è grande.
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
 * Integrazione RK2 (metodo del punto medio) per stabilità numerica:
 *   k1 = f(theta)
 *   k2 = f(theta + 0.5*dt*k1)   ← valuta al punto medio
 *   theta += dt * k2
 *
 * RK2 è più stabile di Eulero per sistemi oscillatori.
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
 * Costruisce il profilo di frequenze naturali omega[KDIM].
 *
 * Profilo base: omega_i = theta0 * decay^i
 *   → decadimento geometrico lungo gli indici (come le catene locali)
 *
 * Bias direzionale tilt = 0.25*(right-left):
 *   omega_i *= (1 + tilt * pos_i)   dove pos_i ∈ [-1, +1]
 *   → right>left: oscillatori con indice alto hanno omega maggiore
 *   → left>right: oscillatori con indice basso hanno omega maggiore
 *
 * right e left influenzano la "direzione" del flusso Kuramoto,
 * non il numero di rotazioni come negli altri opcode.
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
 * Struttura dati dell'opcode
 * ============================================================ */
typedef struct {
  OPDS h;
  MYFLT *out;

  MYFLT *kfreq, *kamp, *itab;
  MYFLT *kstart, *kright, *kleft, *kboth, *ktheta0, *kdecay, *kmove;

  double phase, sr;

  MYFLT base[WT_LEN];     /* wavetable originale */
  MYFLT wt[WT_LEN];       /* wavetable trasformata */

  double theta[KDIM];     /* fasi Kuramoto — PERSISTONO tra k-block!
                           * Questa è la "memoria" del flusso timbrico. */
} OP;

static int load_base(CSOUND *cs, OP *p) {
  FUNC *ft = cs->FTnp2Find(cs, p->itab);
  if (UNLIKELY(ft == NULL || ft->flen < WT_LEN)) return NOTOK;
  for (int i = 0; i < WT_LEN; i++) p->base[i] = ft->ftable[i];
  return OK;
}

/* ============================================================
 * build_wavetable
 * Usa le KDIM fasi Kuramoto correnti come angoli per KDIM rotazioni
 * di Givens non-locali (passo primo) sulla wavetable.
 *
 * i parte da start, avanza di 1 (drift) ad ogni rotazione.
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

/* Loop audio con interpolazione lineare */
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
 * Inizializza fasi Kuramoto a zero: theta_i=0 => wavetable base.
 * ============================================================ */
static int op_init(CSOUND *cs, OP *p) {
  p->phase = 0.0;
  p->sr    = (double)cs->GetSr(cs);
  for (int i = 0; i < KDIM; i++) p->theta[i] = 0.0;
  return OK;
}

/* ============================================================
 * op_perf  [a+k-rate]
 * Ogni k-block:
 *   1) Carica wavetable base
 *   2) Costruisce profilo omega dai parametri
 *   3) Avanza il sistema Kuramoto di dt = 0.002 * kmove
 *   4) Costruisce wavetable con le fasi correnti
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

  /* costruisce profilo omega */
  double omega[KDIM];
  make_omega(omega, theta0, decay, right, left);

  /* avanza il sistema Kuramoto:
   * dt piccolo = evoluzione lenta e stabile
   * kmove=0 => dt=0 => theta fermo => timbro invariato */
  double dt = 0.002 * move;
  kuramoto_step_rk2(p->theta, omega, K, dt);

  /* costruisce la wavetable con le fasi aggiornate */
  build_wavetable(p, start);
  render(p, out, offset, nsmps);
  return OK;
}

/* ============================================================
 * Registrazione opcode
 * Firma: aOut wtgivkura kfreq, kamp, itab,
 *                            kstart, kright, kleft, kboth,
 *                            ktheta0, kdecay, kmove
 * ============================================================ */
static OENTRY localops[] = {
  { "wtgivkura", sizeof(OP), 0, 3, "a", "kkikkkkkkk",
    (SUBR)op_init, (SUBR)op_perf }
};

LINKAGE
