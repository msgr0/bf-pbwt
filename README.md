```bash
make

#for parallel threads number use
OMP_NUM_THREADS=N ./2bfpbwt prs panel.r10.c100000.txt
```

Known issues:
- Buffered version can (will) have issues when input have less column than buffer sizes
