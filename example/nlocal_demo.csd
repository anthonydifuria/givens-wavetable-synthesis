<CsoundSynthesizer>
<CsOptions>
-m0d
-odac
--opcode-dir=./build
--opcode-lib=./build/wtgivnlocal.so
</CsOptions>

<CsInstruments>
sr=48000
ksmps=64
nchnls=2
0dbfs=1

; base: sinusoide
giBase ftgen 1, 0, 1024, 10, 1

instr 1
gkDur = 1
ktrig metro 10
schedkwhen ktrig, 0, 0, 2, 0, gkDur

endin

instr 2
  iHarm = int(rnd(5) + 1)
  kfreq   = 110 * iHarm
  kamp    = 0.1

  kstart  = 0.5

  kright  = 0
  kleft   = 0
  kboth   = 1        ; <-- alza questo per “andare dappertutto”

  ktheta0 = rnd(3.14)        ; più alto = più aggressivo
  kdecay  = 0.9   ; vicino a 1 = mixing lungo
  kmove = 1

  aOut wtgivnlocal kfreq, kamp, 1, kstart, kright, kleft, kboth, ktheta0, kdecay, kmove
  aEnv linseg 0, i(gkDur) / 2, 1, i(gkDur) / 2, 0
  aL, aR pan2 aOut * aEnv, rnd(1)
  outs aL, aR
endin
</CsInstruments>

<CsScore>
i1 0 100
e
</CsScore>
</CsoundSynthesizer>