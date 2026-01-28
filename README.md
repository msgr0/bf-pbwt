```bash
make

#for parallel threads number use
OMP_NUM_THREADS=N ./2bfpbwt prs panel.r10.c100000.txt
```

Known issues:
- Buffered version can (will) have issues when input have less column than buffer sizes

## BM (row X col) roadmap:
- `spr`: needs optimization

## BCF roadmap:
- `lin`: implemented, needs checking
- `bli`: implemented, needs checking
- `bli[s|m]`: not relevant, should not be implemented
- `ars`: implemented, needs checking
- `bar`: not implemented
- `bar[s|m]`: not relevant, should not be implemented
- `prs`: implemented, needs checking
- `bpr`: not relevant, should not be implemented
- `spr`: not implemented, this is very tricky

## ENC (64bit encoding) roadmap:
```bash
xxx=ars
./bcf2enc file.bcf file.em
./2bfpbwt $xxx file.em
```

- `lin`: w.i.p.
- `bli`: ?
- `blis`: not relevant (?)
- `blim`: not relevant (?)
- `ars`: implemented, check last computation of last window
- `bar`: ?
- `bars`: not relevant (?)
- `barm`: not relevant (?)
- `prs`: w.i.p.
- `bpr`: ?
- `spr`: ?
