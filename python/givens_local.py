"""
GIVENS_LOCAL.py
============================================================================
Interactive demo: wtgivlocal — local adjacent chain
on an 8-sample wavetable.

Applies Givens rotations to ADJACENT pairs (k, k+1) starting from kstart.
The angle decays geometrically: theta_t = theta0 * kmove * decay^t
Effect: slow spectral diffusion, localised around kstart.

Three chains (in order):
  both  : bidirectional from start, symmetric alternating expansion
  right : from start, rightward to the end
  left  : pivot on (start, start+1), then goes back toward the left

Sliders:
  start  : chain starting index
  right  : amplitude/density of chain toward the right [0..1]
  left   : amplitude/density of chain toward the left [0..1]
  both   : amplitude/density of bidirectional chain [0..1]
  theta0 : initial angle [radians]
  decay  : geometric decay of angles [0..1]
  kmove  : global scaler (0=base, 1=full)
============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# ──────────────────────────────────────────────
# Core mathematics
# ──────────────────────────────────────────────

def givens(x, i, j, theta):
    """
    In-place Givens rotation in the (i,j) plane.
    Modifies only x[i] and x[j], preserving x[i]^2 + x[j]^2.
    """
    c, s = np.cos(theta), np.sin(theta)
    xi, xj = x[i], x[j]
    x[i] = c*xi - s*xj
    x[j] = s*xi + c*xj

def peak_normalize(x, eps=1e-12):
    """Normalises to absolute peak = 1.0."""
    m = np.max(np.abs(x))
    return x if m < eps else x / m


# ──────────────────────────────────────────────
# Plugin 1: local adjacent chain (wtgivlocal)
# ──────────────────────────────────────────────

def local_chain(x0, start, right, left, both, theta0, decay, kmove=1.0):
    """
    Applies Givens chains on adjacent pairs.
    theta_eff = theta0 * kmove  (kmove=0 => wavetable unchanged)
    """
    x     = x0.astype(float).copy()
    n     = len(x)
    start = int(np.clip(round(start), 0, n-2))
    decay = float(np.clip(decay, 0.0, 1.0))
    theta_eff = float(theta0) * float(np.clip(kmove, 0.0, 1.0))

    def chain_right(theta):
        th = theta
        for k in range(start, n-1):
            givens(x, k, k+1, th)
            th *= decay

    def chain_left(theta):
        # first pair is always the pivot (start, start+1)
        th = theta
        givens(x, start, start+1, th)
        th *= decay
        for k in range(start, 0, -1):
            givens(x, k-1, k, th)
            th *= decay

    def chain_both(theta):
        # t=0: pivot; then alternates left and right
        th = theta
        givens(x, start, start+1, th); th *= decay
        L, R = start-1, start+1
        while L >= 0 or R <= n-2:
            if L >= 0:
                givens(x, L, L+1, th); th *= decay; L -= 1
            if R <= n-2:
                givens(x, R, R+1, th); th *= decay; R += 1

    if both  != 0.0: chain_both (theta_eff * both)
    if right != 0.0: chain_right(theta_eff * right)
    if left  != 0.0: chain_left (theta_eff * left)

    return peak_normalize(x)


# ──────────────────────────────────────────────
# GUI
# ──────────────────────────────────────────────

n   = 8
idx = np.arange(n)
x0  = np.sin(2*np.pi*idx/n)   # base wavetable: sinusoid

p0  = dict(start=3, right=0.0, left=0.0, both=1.0, theta0=0.785, decay=1.0, kmove=1.0)

fig, ax = plt.subplots(1, 1, figsize=(11, 6))
plt.subplots_adjust(left=0.10, right=0.75, top=0.90, bottom=0.36)

ax.set_title("wtgivlocal — local adjacent chain (8 samples)", fontsize=11)
ax.set_ylim(-1.2, 1.2)
ax.set_ylabel("value")
ax.set_xlabel("sample index")
ax.grid(True, alpha=0.3)

y_local = local_chain(x0, **p0)
(base1,) = ax.plot(idx, x0,      marker="o", lw=1, label="base (sinusoid)")
(line1,) = ax.plot(idx, y_local, marker="s", lw=2, label="local transformed")
ax.legend(loc="upper right", fontsize=9)

SL = 0.10
SW = 0.70
def sax(bottom):
    return plt.axes([SL, bottom, SW, 0.028])

ax_kmove  = sax(0.27)
ax_start  = sax(0.23)
ax_right  = sax(0.19)
ax_left   = sax(0.15)
ax_both   = sax(0.11)
ax_theta0 = sax(0.07)
ax_decay  = sax(0.03)

s_kmove  = Slider(ax_kmove,  "kmove",        0.0, 1.0,  valinit=p0["kmove"])
s_start  = Slider(ax_start,  "start",        0,   n-2,  valinit=p0["start"], valstep=1)
s_right  = Slider(ax_right,  "right",        0.0, 1.0,  valinit=p0["right"])
s_left   = Slider(ax_left,   "left",         0.0, 1.0,  valinit=p0["left"])
s_both   = Slider(ax_both,   "both",         0.0, 1.0,  valinit=p0["both"])
s_theta0 = Slider(ax_theta0, "theta0 [rad]", 0.0, 3.14, valinit=p0["theta0"])
s_decay  = Slider(ax_decay,  "decay",        0.0, 1.0,  valinit=p0["decay"])

def update(_):
    params = dict(
        start=int(s_start.val),
        right=float(s_right.val),
        left=float(s_left.val),
        both=float(s_both.val),
        theta0=float(s_theta0.val),
        decay=float(s_decay.val),
        kmove=float(s_kmove.val),
    )
    line1.set_ydata(local_chain(x0, **params))
    fig.canvas.draw_idle()

for s in (s_kmove, s_start, s_right, s_left, s_both, s_theta0, s_decay):
    s.on_changed(update)

plt.show()
