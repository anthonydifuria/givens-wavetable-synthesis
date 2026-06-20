"""
GIVENS_NLOCAL.py
============================================================================
Interactive demo: wtgivnlocal — non-local prime step
on an 8-sample wavetable.

Applies Givens rotations to distant pairs (i, i+PRIME_STEP) % N.
PRIME_STEP=3 (coprime with 8 in the demo, 37 in the real opcode on 1024).
Effect: global spectral mixing from the very first steps.

right/left/both control the NUMBER of rotations (0..max_steps).
Order: non-local (both) -> forward (right) -> backward (left)

Sliders:
  start  : chain starting index
  right  : forward rotation density [0..1]
  left   : backward rotation density [0..1]
  both   : non-local rotation density [0..1]
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
# Plugin 2: non-local prime step (wtgivnlocal)
# ──────────────────────────────────────────────

def nonlocal_chain(x0, start, right, left, both, theta0, decay, kmove=1.0,
                   prime_step=3, max_steps=24):
    """
    Applies Givens rotations to NON-LOCAL pairs.
    prime_step=3 (8-sample demo); in the real opcode on 1024: prime_step=37.
    """
    x     = x0.astype(float).copy()
    n     = len(x)
    start = int(np.clip(round(start), 0, n-1))
    decay = float(np.clip(decay, 0.0, 1.0))
    theta_eff = float(theta0) * float(np.clip(kmove, 0.0, 1.0))

    stepsR = int(np.floor(np.clip(right, 0, 1) * max_steps))
    stepsL = int(np.floor(np.clip(left,  0, 1) * max_steps))
    stepsB = int(np.floor(np.clip(both,  0, 1) * max_steps))

    def apply_right(steps):
        th = theta_eff
        i  = start % n
        for _ in range(steps):
            j = (i + 1) % n
            givens(x, i, j, th)
            th *= decay
            i   = j

    def apply_left(steps):
        th = theta_eff
        i  = start % n
        for _ in range(steps):
            j = (i - 1) % n
            givens(x, i, j, th)
            th *= decay
            i   = j

    def apply_both(steps):
        # prime step: global traversal of the space
        th = theta_eff
        i  = start % n
        for _ in range(steps):
            j = (i + prime_step) % n
            givens(x, i, j, th)
            th *= decay
            i   = (i + 1) % n   # drift

    if stepsB > 0: apply_both (stepsB)
    if stepsR > 0: apply_right(stepsR)
    if stepsL > 0: apply_left (stepsL)

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

ax.set_title("wtgivnlocal — non-local prime step (8 samples)", fontsize=11)
ax.set_ylim(-1.2, 1.2)
ax.set_ylabel("value")
ax.set_xlabel("sample index")
ax.grid(True, alpha=0.3)

y_non = nonlocal_chain(x0, **p0, prime_step=3, max_steps=24)
(base2,) = ax.plot(idx, x0,    marker="o", lw=1, label="base (sinusoid)")
(line2,) = ax.plot(idx, y_non, marker="s", lw=2, label="non-local transformed")
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
    line2.set_ydata(nonlocal_chain(x0, **params, prime_step=3, max_steps=24))
    fig.canvas.draw_idle()

for s in (s_kmove, s_start, s_right, s_left, s_both, s_theta0, s_decay):
    s.on_changed(update)

plt.show()
