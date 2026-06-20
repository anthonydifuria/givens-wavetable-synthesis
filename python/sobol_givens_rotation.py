"""
SOBOL_GIVENS_ROTATION.py
============================================================================
Interactive demo: wavetable navigation with SOBOL 1D-guided angles.

DETERMINISTIC behaviour:
  The Sobol generator is reset (seed=0) at every GUI update,
  exactly as in the C code (sobol_reset() at each k-block).
  → Same parameters  =  same timbre every time
  → kmove=0 returns exactly the base wavetable

What is the Sobol sequence?
  A low-discrepancy quasi-random sequence: points fill [0,1)
  more uniformly than pseudo-random numbers, avoiding clusters and gaps.
  In musical terms: every region of the SO(N) timbral space is
  explored in a balanced way, with no over- or under-represented zones.

  Algorithm (Gray-code increment):
    x_n = x_{n-1} XOR v_c     where c = trailing_zeros(n)
    v_c = 1 << (31 - c)

  The rotation angle at step t is:
    u      ∈ [0,1)    from Sobol
    z      = 2u - 1  ∈ [-1, 1]
    theta_t = kmove * theta0 * decay^t * z

Sliders:
  kmove  : transformation intensity (0=base, 1=full)
  start  : starting index
  right  : forward local mixing density [0..1]
  left   : backward local mixing density [0..1]
  both   : non-local prime-step mixing density [0..1]
  theta0 : angle scale [radians]
  decay  : geometric decay [0..1]
============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


# ──────────────────────────────────────────────
# Deterministic 1D Sobol generator
# ──────────────────────────────────────────────

class Sobol1D:
    """
    1D Sobol with simplified direction numbers: v_c = 1 << (31 - c)
    Identical to the C implementation in the wtgivsobol opcode.

    Usage:
      sob = Sobol1D()          # seed=0
      sob.next()               # → float in [0,1)
      sob.reset()              # → back to the initial point
    """
    def __init__(self, seed=0):
        self.reset(seed)

    def reset(self, seed=0):
        """Deterministic reset: same sequence from the start."""
        self.i = 0
        self.x = np.uint32(seed)

    @staticmethod
    def ctz32(n):
        """Count trailing zero bits of n > 0."""
        c = 0
        while (n & 1) == 0:
            n >>= 1
            c  += 1
        return c

    def next(self):
        """
        Next sample in [0,1).
        Uses the Gray-code increment:
          i += 1
          c  = trailing_zeros(i)
          v  = 1 << (31 - c)
          x ^= v
        """
        self.i += 1
        c  = Sobol1D.ctz32(self.i)
        v  = np.uint32(1 << (31 - c))
        self.x ^= v
        return float(self.x) / float(2**32)


# ──────────────────────────────────────────────
# Core mathematics
# ──────────────────────────────────────────────

def peak_normalize(x, eps=1e-12):
    m = np.max(np.abs(x))
    return x if m < eps else x / m

def givens(x, i, j, theta):
    c, s  = np.cos(theta), np.sin(theta)
    xi, xj = x[i], x[j]
    x[i]  = c*xi - s*xj
    x[j]  = s*xi + c*xj


# ──────────────────────────────────────────────
# Build wavetable with Sobol (DETERMINISTIC)
# ──────────────────────────────────────────────

def build_wavetable_sobol(x0, start, right, left, both, theta0, decay, kmove,
                          prime_step=3, max_steps=24, sobol_seed=0):
    """
    Builds the transformed wavetable with Sobol-guided angles.

    RESET at the start: same parameters → same wavetable (as in the C opcode).

    Application order (identical to the opcode):
      1. both  : non-local prime-step mixing
      2. right : forward local mixing
      3. left  : backward local mixing

    Angle at step t:
      u = Sobol.next()  ∈ [0,1)
      z = 2u - 1        ∈ [-1,1]
      theta_t = kmove * theta0 * decay^t * z
    """
    n      = len(x0)
    x      = x0.astype(float).copy()
    start  = int(np.clip(round(start), 0, n-1))
    right  = float(np.clip(right, 0, 1))
    left   = float(np.clip(left,  0, 1))
    both   = float(np.clip(both,  0, 1))
    theta0 = float(theta0)
    decay  = float(np.clip(decay, 0, 1))
    move   = float(np.clip(kmove, 0, 1))

    stepsR = int(np.floor(right * max_steps))
    stepsL = int(np.floor(left  * max_steps))
    stepsB = int(np.floor(both  * max_steps))

    # Deterministic RESET (like sobol_reset() in C)
    sob = Sobol1D(seed=sobol_seed)

    def theta(t):
        u = sob.next()           # [0,1)
        z = 2.0*u - 1.0          # [-1,1] centred on zero
        return move * theta0 * (decay**t) * z

    t = 0   # global step counter

    # 1) both: non-local prime step
    i = start % n
    for _ in range(stepsB):
        j = (i + prime_step) % n
        givens(x, i, j, theta(t)); t += 1
        i = (i + 1) % n          # drift

    # 2) right: forward local
    i = start % n
    for _ in range(stepsR):
        j = (i + 1) % n
        givens(x, i, j, theta(t)); t += 1
        i = j

    # 3) left: backward local
    i = start % n
    for _ in range(stepsL):
        j = (i - 1) % n
        givens(x, i, j, theta(t)); t += 1
        i = j

    return peak_normalize(x)


# ──────────────────────────────────────────────
# GUI
# ──────────────────────────────────────────────

n   = 8
idx = np.arange(n)
x0  = np.sin(2*np.pi*idx/n)

p0  = dict(start=3, right=0.0, left=0.0, both=1.0, theta0=0.785, decay=1.0, kmove=1.0)

fig, ax = plt.subplots(1, 1, figsize=(11, 6))
plt.subplots_adjust(left=0.10, right=0.75, top=0.90, bottom=0.36)

ax.set_title("wtgivsobol — Sobol-driven Givens (8 samples)\n"
             "kmove=0 → base sinusoid unchanged", fontsize=11)
ax.set_ylim(-1.2, 1.2)
ax.set_xlabel("sample index")
ax.set_ylabel("value")
ax.grid(True, alpha=0.3)

y = build_wavetable_sobol(x0, **p0, prime_step=3, max_steps=24, sobol_seed=0)
(base_line,) = ax.plot(idx, x0, marker="o", lw=1, label="base")
(line,)      = ax.plot(idx, y,  marker="s", lw=2, label="sobol transformed")
ax.legend(loc="upper right", fontsize=9)

SL, SW = 0.10, 0.62
def sax(b): return plt.axes([SL, b, SW, 0.028])

ax_kmove  = sax(0.30)
ax_start  = sax(0.26)
ax_right  = sax(0.22)
ax_left   = sax(0.18)
ax_both   = sax(0.14)
ax_th0    = sax(0.10)
ax_decay  = sax(0.06)

s_kmove  = Slider(ax_kmove, "kmove",         0.0, 1.0,  valinit=p0["kmove"])
s_start  = Slider(ax_start, "start",         0,   n-1,  valinit=p0["start"], valstep=1)
s_right  = Slider(ax_right, "right",         0.0, 1.0,  valinit=p0["right"])
s_left   = Slider(ax_left,  "left",          0.0, 1.0,  valinit=p0["left"])
s_both   = Slider(ax_both,  "both",          0.0, 1.0,  valinit=p0["both"])
s_th0    = Slider(ax_th0,   "theta0 [rad]",  0.0, 8.0,  valinit=p0["theta0"])
s_decay  = Slider(ax_decay, "decay",         0.0, 1.0,  valinit=p0["decay"])

def update(_):
    params = dict(
        start=int(s_start.val),
        right=float(s_right.val),
        left=float(s_left.val),
        both=float(s_both.val),
        theta0=float(s_th0.val),
        decay=float(s_decay.val),
        kmove=float(s_kmove.val),
    )
    y = build_wavetable_sobol(x0, **params, prime_step=3, max_steps=24, sobol_seed=0)
    line.set_ydata(y)
    fig.canvas.draw_idle()

for s in (s_kmove, s_start, s_right, s_left, s_both, s_th0, s_decay):
    s.on_changed(update)

plt.show()
