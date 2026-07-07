# givens-wavetable-synthesis

Wavetable navigation in high-dimensional space via Givens rotations.

A discrete audio signal of N samples is treated as a point in R^N. Givens
rotations act on pairs of coordinates through the sine and cosine of an angle,
displacing the point along the surface of a hypersphere while preserving the
L² norm. Four strategies for selecting rotation pairs and angles are provided,
each exploring different regions of the waveform space.

---

## Strategies

| Opcode / Script | Strategy | Pair selection | Angle sequence |
|---|---|---|---|
| `wtgivlocal` / `givens_local.py` | Local adjacent | adjacent (k, k+1) | geometric decay from pivot |
| `wtgivnlocal` / `givens_nlocal.py` | Non-local prime step | (i, (i+p) mod N), p=37 | geometric decay |
| `wtgivsobol` / `sobol_givens_rotation.py` | Sobol quasi-random | non-local prime step | low-discrepancy Sobol sequence |
| `wtgivkura` / `kuramoto_givens_rotation.py` | Kuramoto coupled oscillators | non-local prime step | phases of a Kuramoto ODE system |

All strategies preserve signal energy (L² norm). Peak normalization is applied
after each rotation sequence.

---

## Repository structure

```
givens-wavetable-synthesis/
├── csound/
│   ├── wtgivlocal/      # local adjacent chain opcode
│   ├── wtgivnlocal/     # non-local prime-step opcode
│   ├── wtgivsobol/      # Sobol quasi-random opcode
│   └── wtgivkura/       # Kuramoto coupled-oscillator opcode
├── python/              # interactive prototypes (N=8, matplotlib GUI)
├── examples/            # Csound .csd usage examples
├── .gitignore
└── README.md
```

---

## Csound opcodes

### Requirements

- Csound ≥ 6.18
- CMake ≥ 3.15
- C compiler (gcc or clang)

### Build (one opcode)

```bash
cd csound/wtgivlocal
mkdir build && cd build
cmake ..
make
```

The compiled `.so` is placed in `build/`. Copy it to your Csound plugin directory
or pass `-+plugins_dir=.` when running Csound.

### Build all opcodes

```bash
for op in wtgivlocal wtgivnlocal wtgivsobol wtgivkura; do
  cd csound/$op
  mkdir -p build && cd build
  cmake .. && make
  cd ../../..
done
```

### Opcode syntax (common to all four)

```csound
aOut  opcodename  kfreq, kamp, itab,
                  kstart, kright, kleft, kboth,
                  ktheta0, kdecay, kmove
```

| Parameter | Local / NLocal / Sobol | Kuramoto |
|---|---|---|
| `kstart` | pivot index | pivot index |
| `kright` | right chain density [0–1] | ω-profile bias R [0–1] |
| `kleft` | left chain density [0–1] | ω-profile bias L [0–1] |
| `kboth` | bidirectional density [0–1] | coupling strength K [0–10] |
| `ktheta0` | initial angle [rad] | ω scale [rad] |
| `kdecay` | geometric decay [0–1] | ω decay [0–1] |
| `kmove` | intensity [0–1] | evolution speed [0–10] |

`kmove = 0` always returns the unmodified wavetable.

### Minimal example

```csound
giSin  ftgen  0, 0, 1024, 10, 1

instr 1
  kfreq  =  220
  kamp   =  0.5
  aOut   wtgivlocal  kfreq, kamp, giSin, \
                     3, 0.0, 0.0, 1.0,   \
                     0.785, 0.9, line:k(0, p3, 1)
  outs   aOut, aOut
endin
```

`kmove` evolves from 0 to 1 over the note duration, producing a gradual
transformation from the original sinusoid to the rotated wavetable.

---

## Python prototypes

Each script runs an interactive matplotlib GUI operating on N=8 samples
with a base sinusoid x_n = sin(2πn/8).

### Requirements

```bash
pip install numpy matplotlib
```

### Run

```bash
python python/givens_local.py
python python/givens_nlocal.py
python python/sobol_givens_rotation.py
python python/kuramoto_givens_rotation.py
```

Parameters are controlled via sliders. For the Kuramoto prototype, the `[STEP]`
button advances the system by one block; `[RESET]` sets all phases to zero.

---

## Mathematical background

A wavetable of N samples is a vector **x** ∈ R^N. Orthogonal transformations
that preserve ‖**x**‖₂ form the rotation group SO(N). A Givens rotation
G(i, j, θ) acts on two components only:

```
x'_i =  cos θ · x_i − sin θ · x_j
x'_j =  sin θ · x_i + cos θ · x_j
```

Any rotation in R^N can be expressed as a composition of Givens rotations.
The choice of which pairs (i, j) to rotate and which angles θ to apply
determines the trajectory through the waveform space.

For N = 1024 the rotation group has dimension N(N−1)/2 = 523,776 degrees
of freedom.

---

## License

MIT

---

## Authors

Luca Bimbi

Anthony Di Furia — LEAP, Laboratorio ElettroAcustico Permanente, Rome - Conservatorio N. Piccinni, Bari 

Giuseppe Silvi — Conservatorio A. Casella, L'Aquila
