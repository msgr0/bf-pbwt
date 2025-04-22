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
- `bpr`: not implemented
- `spr`: not implemented, this is very tricky

## BM (modified preprocess) roadmap:
- `lin`: ?
- `bli`: ?
- `blis`: ?
- `blim`: ?
- `ars`: ?
- `bar`: ?
- `bars`: ?
- `barm`: ?
- `prs`: ?
- `bpr`: ?
- `spr`: ?
