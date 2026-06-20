"""
WTGIVKURA.py
============================================================================
Interactive demo: wtgivkura — Kuramoto-driven Givens wavetable navigator.
Mirrors wtgivkura.c exactly:

  - KDIM = 64 Kuramoto oscillators / Givens rotations
  - PRIME_STEP = 37 (coprime with 1024; scaled to 3 for 8-sample demo)
  - dt = 0.002 * kmove  (kmove=0 → system frozen → static timbre)
  - Integration: RK2 (midpoint)
  - omega profile: base_i = theta0 * decay^i, tilt = 0.25*(right-left)*pos_i
  - Givens schedule: non-local pair (i, i+PRIME_STEP), drift i+=1 each step
  - Phases theta[] persist across k-blocks → continuous timbral trajectory

Fundamental difference vs Sobol:
  Sobol    → SPATIAL:  fixed params = fixed timbre (static map)
  Kuramoto → TEMPORAL: fixed params = continuous TRAJECTORY (flow)

Syntax (mirrors Csound):
  wtgivkura kfreq, kamp, itab,
            kstart, kright, kleft, kboth,
            ktheta0, kdecay, kmove

Controls:
  [STEP]  : advances Kuramoto by one k-block (dt = 0.002 * kmove)
  [RESET] : resets all phases to zero → base wavetable
  start   : starting index in the wavetable
  right   : omega bias toward the right [0..1]
  left    : omega bias toward the left  [0..1]
  both (K): Kuramoto coupling strength  [0..10]
  theta0  : scale of natural frequencies omega_i
  decay   : geometric decay of omega profile [0..1]
  kmove   : evolution speed — dt = 0.002 * kmove  [0..10]
============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# ──────────────────────────────────────────────
# Constants (matching wtgivkura.c)
# ──────────────────────────────────────────────

N          = 8       # wavetable length (demo; C uses 1024)
KDIM       = 8       # Kuramoto oscillators (demo; C uses 64)
PRIME_STEP = 3       # coprime with N=8 (C uses 37 for N=1024)

# ──────────────────────────────────────────────
# Core mathematics
# ──────────────────────────────────────────────

def givens_inplace(x, i, j, theta):
    """x'[i] = cos*x[i] - sin*x[j],  x'[j] = sin*x[i] + cos*x[j]"""
    c, s = np.cos(theta), np.sin(theta)
    xi, xj = x[i], x[j]
    x[i] = c*xi - s*xj
    x[j] = s*xi + c*xj

def peak_normalize(x, eps=1e-12):
    m = np.max(np.abs(x))
    return x if m < eps else x / m

# ──────────────────────────────────────────────
# Kuramoto model (mirrors kuramoto_rhs + kuramoto_step_rk2 in .c)
# ──────────────────────────────────────────────

def kuramoto_rhs(theta, omega, K):
    """
    dtheta_i = omega_i + (K/N) * Σ_j sin(theta_j - theta_i)
    """
    diff     = theta[None, :] - theta[:, None]   # diff[i,j] = theta_j - theta_i
    coupling = np.sum(np.sin(diff), axis=1) / KDIM
    return omega + K * coupling

def kuramoto_step_rk2(theta, omega, K, dt):
    """
    RK2 midpoint integration (matches wtgivkura.c):
      k1  = f(theta)
      k2  = f(theta + 0.5*dt*k1)
      theta += dt * k2
    """
    k1  = kuramoto_rhs(theta, omega, K)
    mid = theta + 0.5 * dt * k1
    k2  = kuramoto_rhs(mid, omega, K)
    return theta + dt * k2

# ──────────────────────────────────────────────
# Omega profile (mirrors make_omega in .c)
# ──────────────────────────────────────────────

def make_omega(theta0, decay, right, left):
    """
    omega_i = theta0 * decay^i * (1 + tilt * pos_i)
    tilt = 0.25 * (right - left),  pos_i ∈ [-1, +1]
    """
    tilt = 0.25 * (right - left)
    base = np.array([theta0 * (decay**i) for i in range(KDIM)])
    pos  = np.linspace(-1.0, 1.0, KDIM)
    return base * (1.0 + tilt * pos)

# ──────────────────────────────────────────────
# Wavetable builder (mirrors build_wavetable in .c)
# ──────────────────────────────────────────────

def build_wavetable(base, theta, start):
    """
    Apply KDIM non-local Givens rotations using current Kuramoto phases.
    Pair: (i, i+PRIME_STEP) % N, drift i += 1 each step.
    """
    x = base.astype(float).copy()
    i = int(start) % N
    for k in range(KDIM):
        j = (i + PRIME_STEP) % N
        givens_inplace(x, i, j, theta[k])
        i = (i + 1) % N
    return peak_normalize(x)

# ──────────────────────────────────────────────
# GUI
# ──────────────────────────────────────────────

idx  = np.arange(N)
base = np.sin(2*np.pi*idx/N)   # base wavetable: sinusoid

# Kuramoto state — persists across steps (mirrors p->theta[] in .c)
theta_state = np.zeros(KDIM, dtype=float)

fig, ax = plt.subplots(1, 1, figsize=(11, 6))
plt.subplots_adjust(left=0.10, right=0.75, top=0.90, bottom=0.38)

ax.set_title("wtgivkura — Kuramoto flow → Givens angles → wavetable (8 samples)\n"
             "Press [STEP] to advance one k-block  |  kmove=0 → frozen timbre",
             fontsize=11)
ax.set_ylim(-1.2, 1.2)
ax.set_xlabel("sample index")
ax.set_ylabel("value")
ax.grid(True, alpha=0.3)

(base_line,) = ax.plot(idx, base, marker="o", lw=1, label="base")
y = build_wavetable(base, theta_state, start=0)
(line,) = ax.plot(idx, y, marker="s", lw=2, label="kuramoto-wavetable")
ax.legend(loc="upper right", fontsize=9)

# Sliders
SL, SW = 0.10, 0.62
def sax(b): return plt.axes([SL, b, SW, 0.028])

ax_start  = sax(0.31)
ax_right  = sax(0.27)
ax_left   = sax(0.23)
ax_both   = sax(0.19)
ax_th0    = sax(0.15)
ax_decay  = sax(0.11)
ax_kmove  = sax(0.07)

s_start  = Slider(ax_start, "start",            0,    N-1,  valinit=0,   valstep=1)
s_right  = Slider(ax_right, "right (ω bias)",   0.0,  1.0,  valinit=0.2)
s_left   = Slider(ax_left,  "left  (ω bias)",   0.0,  1.0,  valinit=0.2)
s_both   = Slider(ax_both,  "both (K)",          0.0, 10.0,  valinit=1.2)
s_th0    = Slider(ax_th0,   "theta0 (ω scale)",  0.0, 10.0,  valinit=3.0)
s_decay  = Slider(ax_decay, "decay",             0.0,  1.0,  valinit=0.95)
s_kmove  = Slider(ax_kmove, "kmove  [dt=0.002·kmove]", 0.0, 10.0, valinit=1.0)

# Buttons
ax_btn  = plt.axes([0.10, 0.02, 0.20, 0.04])
ax_rst  = plt.axes([0.32, 0.02, 0.20, 0.04])
btn     = Button(ax_btn, "STEP",      hovercolor='0.85')
btn_rst = Button(ax_rst, "RESET θ=0", hovercolor='0.85')

def get_omega():
    return make_omega(
        theta0 = float(s_th0.val),
        decay  = float(s_decay.val),
        right  = float(s_right.val),
        left   = float(s_left.val),
    )

def do_step(_=None):
    """
    Advances Kuramoto by one k-block.
    dt = 0.002 * kmove  — identical to wtgivkura.c op_perf().
    kmove=0 → dt=0 → theta unchanged → static timbre.
    """
    global theta_state
    K     = float(np.clip(s_both.val,  0.0, 10.0))
    move  = float(np.clip(s_kmove.val, 0.0, 10.0))
    dt    = 0.002 * move          # ← exact match with .c
    omega = get_omega()

    theta_state = kuramoto_step_rk2(theta_state, omega, K, dt)

    start = int(s_start.val)
    y = build_wavetable(base, theta_state, start)
    line.set_ydata(y)
    fig.canvas.draw_idle()

def do_reset(_=None):
    """Resets all Kuramoto phases to zero → base wavetable (mirrors op_init)."""
    global theta_state
    theta_state = np.zeros(KDIM, dtype=float)
    y = build_wavetable(base, theta_state, int(s_start.val))
    line.set_ydata(y)
    fig.canvas.draw_idle()

btn.on_clicked(do_step)
btn_rst.on_clicked(do_reset)

plt.show()