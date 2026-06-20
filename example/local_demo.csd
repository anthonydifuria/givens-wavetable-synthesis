<CsoundSynthesizer>
<CsOptions>
-odac
--opcode-dir=./build
--opcode-lib=./build/wtgivlocal.so
</CsOptions>

<CsInstruments>
sr=48000
ksmps=64
nchnls=2
0dbfs=1

; wavetable base 1024 (sine)
giBase ftgen 1, 0, 1024, 10, 1

instr 1
  kfreq   = 300
  kamp    = 0.2

  kstart  = 0.5

  kright  = 0.0
  kleft   = 0.0
  kboth   = 1.0

  kT  linseg 0, 10, 1
  ktheta0 =  3.14 * 1
  ;printk 0.1, ktheta0
  kdecay  = 0.9        ; scaling/decay (vicino a 1 = dispersione lunga)
  kmove = 1.0

  aOut wtgivlocal kfreq, kamp, 1, kstart, kright, kleft, kboth, ktheta0, kdecay, kmove
  outs aOut, aOut
endin
</CsInstruments>

<CsScore>
i1 0 100
e
</CsScore>
</CsoundSynthesizer>