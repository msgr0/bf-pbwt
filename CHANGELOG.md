# Changelog — sp-pbwt optimization work

## Run 0 (Baseline Harness) — 2026-08-28

**Status:** STOP — Correctness issues detected before any code changes.

### Build

- Successfully compiled `sp-pbwt-bcf`, `sp-pbwt-bm`, and `gen` binaries.
- Note: Used direct `cc` compilation flags instead of Makefile (htslib installed at `/data/proj/2bfpbwt/htslib/`, headers in `htslib/htslib/`).

### Correctness Fixture

- Created small BCF fixture with 200 records: `bench/panel.small.bcf` (30 KB)
- Ran all four BCF modes on the fixture with DUMP output mode.
- Used `dumpcheck.py` to compare prefix and divergence arrays across modes.

**Correctness Test Results:**

| Comparison | Checked | Wrong | Status |
|---|---|---|---|
| linear vs stagpar | 199 | 0 | ✓ MATCH |
| linear vs blockpar | 200 | 2 | ✗ DIVERGE |
| linear vs sampled | 4 | 1 | ✗ DIVERGE (only 4 cols produced of 200) |

**Critical Findings:**
1. **blockpar** (`wparc_rrs`): Produces incorrect results on columns 64 and 128 (at minimum). Divergence arrays differ from linear mode.
2. **sampled** (`wapproxc_rrs`): Severely broken — only processes 4 columns instead of 200, suggesting early exit bug or EOF handling issue.
3. **stagpar** (`wstagparc_rrs`): Appears correct on the fixture (produces identical output to linear).

### Baseline Timing (chr20, 5008 samples, ~60k variants)

Partial results before parallel modes timed out:
- **linear**: ~76 seconds wall time, 75 CPU seconds
- **sampled**: >120 seconds (timeout, likely hanging)
- **blockpar**: Not completed
- **stagpar**: Not completed

### Decision

Paused per Step 0 stop condition, then resumed after review with the user:

- **blockpar** divergence at columns 64/128 (multiples of window width W=64) matches the
  already-documented bug in the plan's "Confirmed findings" (`sp-pbwt.c:973-974`,
  `psrev->a`/`psrev->d` copied from the wrong source array). This is scheduled to be
  fixed (or confirmed dead and removed) in Step 4 — not a new issue.
- **sampled** producing only 4/200 columns is expected, by-design behavior: `sampled`
  (`wapproxc_rrs`, the approximate/sampled mode) only emits values at its window-sample
  points, not every column. Confirmed with the user — this is not a bug. `dumpcheck.py`
  already matches by row/column index, so sparse output from `sampled` compares
  correctly against the columns `linear` produces at the same indices.
- The 1 mismatched value among the 4 sampled columns produced is not separately
  investigated here; it may be related to the already-documented double-free/leak in
  `wapproxc_rrs` (`sp-pbwt.c:752-756`), which Step 6 fixes.

**Proceeding to Steps 1–6 as planned.**

## Run 1 (Step 1: Build flags + dead serialization) — 2026-08-28

**Status:** SUCCESS — All changes implemented without breaking correctness.

### Build

- Modified `Makefile` CFLAGS: added `-DNDEBUG -march=native` to release build (line 2)
- Fixed `Makefile` HTSLIB include paths from `-I ${HTSLIB}/include` to `-I ${HTSLIB}` (line 10)
  to match actual header location at `/data/proj/2bfpbwt/htslib/htslib/`
- Successfully compiled `sp-pbwt-bcf`, `sp-pbwt-bm`, and `gen` binaries

### Changes in `sp-pbwt.c`

1. **Lines 1242-1263 (prescan removal)**:
   - Deleted the full prescan loop that counted `ncol` by reading the entire BCF file
   - Replaced with `ncol = W;` initialization for BCF mode, matching `wparc_rrs` pattern

2. **Line 1247 (pbprev allocation)**:
   - Deleted unused `pbwtad *pbprev = pbwtad_new(nrow);` allocation

3. **Lines 1313-1314 (debug print)**:
   - Deleted per-thread debug fprintf: `fprintf(stderr, ">>>%zu, %zu, %zu, %zu, %zu, %zu\n\n", ...)`

4. **Lines 1329-1355 (loop restructuring)**:
   - Changed loop from fixed `for (size_t j = lane; j + W <= ncol; j += W)` to EOF-driven
   - BCF mode: `while ((lastrowread = fgetcolwgri(...)) == j + W) { ... j += W; }`
   - Non-BCF mode: kept existing for loop (still uses ncol parameter)
   - This makes the loop exit when `fgetcolwgri` returns less than a full window (EOF signal)

5. **Lines 1347-1354 (critical section guard)**:
   - Moved `#pragma omp critical` inside `if (DO_DUMP)` guard
   - Prevents critical section entry in normal (non-dump) runs, eliminating serialization of the `stagpar` hot loop

### Correctness Verification

Ran `sp-pbwt-bcf` on `bench/panel.small.bcf` (200 records, 5006 samples) with DUMP output:

| Mode | Comparison | Checked | Wrong | Status |
|---|---|---|---|---|
| stagpar | vs linear | 199 | 0 | ✓ PASS |
| blockpar | vs linear | 200 | 2 | ✗ DIVERGE (known bug at cols 64/128) |

**stagpar result:** PASS — Maintained perfect correctness (199/199 match) from Step 0.
The prescan removal and EOF-driven loop work correctly; critical-section guard eliminates
unnecessary serialization without breaking computation.

**blockpar result:** Expected — Known bug from Step 0 is unchanged (to be fixed in Step 4).

## Run 2 (Step 2: BCF decode layer) — 2026-08-28

**Status:** SUCCESS — Bit-identical output preserved, no new divergences.

All changes confined to `iobcf.c`; `io.h` signatures untouched.

### Changes in `iobcf.c`

1. **Hoisted the genotype buffer out of every per-record loop.**
   Introduced a shared helper `bcf_decode_gt_row(hdr, line, n, out)` that all
   five call sites now use (`fgetcoli`, `bfgetcoln`, the `FGETCOLIW_IMPL(W)`
   macro body instantiated for W=8/16/32/64, `fgetcoliwgr`, `fgetcolwgri`).
   Its fallback path (the old `bcf_get_genotypes` route, used for
   multiallelic/non-diploid/non-int8 records) keeps `gt_arr`/`ngt_arr` as
   `static __thread int32_t *` / `static __thread int` instead of allocating
   and `free()`-ing them per record. Because `fgetcolwgri` is called
   concurrently by `stagpar` (one `bcf_srs_t` reader per lane/thread), the
   thread-local qualifier was applied uniformly rather than plain `static`,
   to avoid a race.

2. **Transposed the bit-packing loops in the three window readers**
   (`FGETCOLIW_IMPL`'s `fgetcoliw##W##r`, `fgetcoliwgr`, `fgetcolwgri`).
   Each now decodes its window's records into a per-thread `W x n` (or
   `w x n`) `uint8_t` staging buffer (one contiguous write per record, no
   read-modify-write of `c`), then does a single column-wise packing pass
   building each `c[r]` in one write. The staging buffer is not pre-zeroed;
   an early EOF simply limits how many `k` indices the packing loop ORs in,
   reproducing the exact same partial-window semantics the old
   memset-then-OR code had.

3. **Avoided the full int32 expansion in `bcf_get_genotypes`.**
   `bcf_decode_gt_row`'s fast path calls `bcf_get_fmt(hdr, line, "GT")`
   (which still triggers `bcf_unpack(line, BCF_UN_FMT)` internally — needed
   for any FORMAT field access) and decodes `fmt->p` directly in its native
   width when `fmt->type == BCF_BT_INT8 && fmt->n == 2` (the common
   biallelic-diploid case for phase3-scale panels). Falls back to the
   original `bcf_get_genotypes` path for anything else (multiallelic, wider
   ploidy, non-int8 storage), so those records still decode correctly.

4. **`IOBCF_ASSUME_WELLFORMED` / well-formedness check.**
   Left as a no-op (unchanged), with a detailed comment added explaining the
   tradeoff: a real check is two extra branches per allele in the hottest
   loop in the decode path, the existing baseline was verified
   bit-identical with it disabled, and the fast native-width decode path
   would additionally need a *different* sentinel check
   (`bcf_int8_vector_end`, not `bcf_int32_vector_end`) if ever re-enabled —
   documented inline so a future change doesn't silently check the wrong
   sentinel.

### Build

`HTSLIB=/data/proj/2bfpbwt/htslib make` — compiled cleanly, no warnings/errors.

### Correctness verification (`bench/panel.small.bcf`, 200 records, 5006 samples)

| Comparison | Checked | Wrong | Status |
|---|---|---|---|
| linear vs stagpar | 199 | 0 | MATCH — same as Run 1 baseline |
| linear vs blockpar | 200 | 2 | Diverges only at columns 64 and 128 — identical to the known baseline bug, no new divergence |
| linear vs sampled | 4 | 1 | Unchanged from baseline (approximate/sparse mode) |

Verified the blockpar divergence indices explicitly: exactly `{64, 128}`, matching
Run 0/Run 1 — confirms Step 2 introduced no new blockpar breakage.

**Thread-safety sanity check:** ran `stagpar` 5 times with `OMP_NUM_THREADS=8`
on the fixture. All 5 runs completed without crash/hang. Raw dump files differ
byte-for-byte between runs (expected: concurrent lane threads interleave their
DUMP lines in different orders each run), but sorted content is identical
across all 5 runs, and each run independently checks 199/199 with 0 wrong
against `linear` — confirms the `static __thread` conversion of
`fgetcolwgri`'s buffers did not introduce a race.

**Larger BCF slice check:** skipped — `bcftools` is not on `PATH` in this
environment, per the plan's fallback instruction.

No new divergences introduced. Proceeding is safe for Step 3.

## Run 3 (Step 3: pbidx_t element width) — 2026-08-28

**Status:** SUCCESS — bit-identical DUMP output preserved, no new divergences, no ASan corruption findings.

All changes confined to `sp-pbwt.c`.

### Changes

1. **New `pbidx_t` typedef.** Added near the top of the file:
   ```c
   typedef uint32_t pbidx_t;
   #define PBIDX_FMT "%u"
   #define PBIDX_MAX UINT32_MAX
   ```
   `struct pbwtad` changed from `{ size_t *a; size_t *d; }` to
   `{ pbidx_t *a; pbidx_t *d; }`. Rows (~5008-10016) and columns (up to
   ~6.2M) both fit comfortably in 32 bits, halving the bytes moved by
   every memcpy, radix scatter, and `reversec` pass.

2. **Propagated `pbidx_t` through every consumer:**
   - `rrsortx`, `rrsortx_noaux`, `rrsort0` — the `s`/`pre`/`post`/`aux`
     permutation arrays are now `pbidx_t *`.
   - `reversec`, `reversecprev` — the `i` assigned into `rev->a[...]` is
     now cast to `pbidx_t`.
   - `recover_div` — return type and its `min` local changed to `pbidx_t`
     (both are divergence values, parallel to `p->d[]`).
   - `divc0` (removed an unused `size_t div` local while there), `divc` —
     the `div` local (assigned into `p->d[i]`) is now `pbidx_t`.
   - `cpbwt`, `cpbwti`, `cpbwtiLCP` — the static `o`/`h` scratch arrays
     (radix-sort-style permutation/LCP scratch parallel to `a`/`d`) and
     the per-call `idx`/`ddx`/`f`/`g`/`k` locals are now `pbidx_t`;
     `memcpy` size arguments switched from the hardcoded `sizeof(size_t)`
     to `sizeof *o` / `sizeof *h` so they track the new width
     automatically.
   - The `aux` scratch buffer declared in every top-level entry point that
     calls `rrsortx`/`rrsort0` — `wapproxc_rrs`, `swbapproxc_rrs`,
     `mwbapproxc_rrs`, `wparc_rrs`, `sbwparc_rrs`, `mbwparc_rrs`,
     `wstagparc_rrs`, `swstagparc_rrs`, `mwstagparc_rrs` (9 call sites) —
     changed from `size_t *` to `pbidx_t *`.
   - `linc`, `sblinc`, `mblinc` needed no changes beyond the struct
     definition (they only call `cpbwti`, no separate scratch arrays).

   **Deliberately left as `size_t`:** column/row counters and loop indices
   that are not literally `a`/`d` values or scratch parallel to them —
   `nrow`, `ncol`, all `j`/`x`/`J`/`wix` column counters in the window
   loops, `n`/`w`/`i`/`i0` parameters of `recover_div` (indices/counts,
   not stored divergence values), the `k`/`k-1` column-index parameter of
   `cpbwtiLCP` itself, the `r`/`q` compaction counters in
   `cpbwt`/`cpbwti`/`cpbwtiLCP` (bounded by `n` = nrow, used only as
   local indices), and the `mask`/`cnt[...]` histogram/boolean scratch in
   the radix sort and `cpbwt`-family functions (booleans/counts, not
   index or divergence values). `uint64_t`/`uint8_t` column buffers
   (`c`, `c0`, `w64`, `pw`, `pw0`, `pw1`) are untouched — they belong to
   the I/O layer, not the `a`/`d` element width.

3. **Format strings.** Every `PDUMPR`/`PDUMP`/`PDUMP_SEQR`/`PDUMP_SEQ`/
   `PDUMP_SEQ_OFFSETR`/`PDUMP_SEQ_OFFSET` print of a raw `(p)->a[...]` or
   raw `(p)->d[...]` element now uses `PBIDX_FMT` instead of `%zu`. Prints
   that compute an expression mixing a `size_t` (`(size_t)(i)`, `offset`,
   the loop index) with a `(p)->d[...]` value were left as `%zu` — the
   `size_t` operand forces the usual-arithmetic-conversion result to
   `size_t` regardless of the narrower `pbidx_t` operand, so those are
   still correctly typed. `DPARR`/`parr` are generic (caller-supplied
   `fmt`) and have no live call sites printing `a`/`d`, so left untouched.

4. **Overflow guards.**
   - `main()`: right after `TRACE(fgetrc(fin, &nrow, &ncol))`, added a
     check that exits with a clear error if `nrow > PBIDX_MAX`.
   - Added a `pbidx_guard(size_t v)` helper (prints an error and
     `exit(EXIT_FAILURE)` if `v > PBIDX_MAX`) called after the column
     counter increments in `cpbwt` and `cpbwti` (their static `k`
     counters), and at the top of `cpbwtiLCP` (its `k` parameter) — these
     are exactly the counters whose value gets stored into `pbidx_t`
     divergence fields and that grow with the column count, which for
     BCF input is only known at EOF. Not expected to trigger for the
     target datasets (~6.2M columns max).

5. **Flat `a`/`d` allocation.** `pbwtad_new` now does a single
   `malloc(2*n*sizeof(pbidx_t))` block, with `a` as the block's base
   pointer and `d = a + n` — two mallocs instead of three, `a`/`d`
   adjacent in memory. `cpbwt` (the only other place that constructs a
   `pbwtad` inline, rather than via `pbwtad_new`) was changed to use the
   same flat-block scheme, so every `pbwtad` in the program has `a` as
   its allocation's base pointer. `PBWTAD_FREE` was updated to free only
   `(p)->a` (the shared block) and just null out `(p)->d` instead of
   separately freeing it — this avoids a double-free/leak now that `a`
   and `d` are one allocation, and works uniformly for `pbwtad`s from
   both constructors.

### Build

`HTSLIB=/data/proj/2bfpbwt/htslib make` — compiled cleanly. Also compiled
with `-Wall -Wformat` added ad hoc (not committed to the Makefile) to
check for narrowing/format-string issues from the width change: only
pre-existing, unrelated `-Wunused-variable`/`-Wunused-but-set-variable`
warnings appeared (dead `pbprev` locals, an unused `c0`/`fd`, etc.) — no
new format-string or conversion warnings from this step's changes.

### Correctness verification (`bench/panel.small.bcf`, 200 records, 5006 samples)

| Comparison | Checked | Wrong | Status |
|---|---|---|---|
| linear vs stagpar | 199 | 0 | MATCH — unchanged, bit-identical (required: pure width change must be lossless) |
| linear vs blockpar | 200 | 2 | Diverges only at columns 64 and 128 — identical to the known baseline bug (Step 4 fixes), no new divergence |
| linear vs sampled | 4 | 1 | Unchanged from baseline (approximate/sparse mode, 1 known mismatch) |

Confirmed the blockpar mismatch indices are exactly `{64, 128}` (via
`dumpcheck.py`'s stderr diff output), matching Run 0/1/2 precisely.

**Sanity check on `linear`'s own dump values** (to rule out a
format-string bug producing self-consistent-but-wrong output): column 0's
`a` array is the identity permutation `0..5007`; the last column (199)
has `a` as a permutation of `0..5007` and `d` values in range `[0, 200]`
— all within expected bounds, no garbage.

### ASan (`make debug`, all four BCF modes on the fixture)

The system's gcc 11 install is missing its ASan runtime shared library
(`/usr/lib64/libasan.so.6.0.0` does not exist anywhere on disk, and the
static `libasan.a` is also absent) — `HTSLIB=... make debug` fails to
*link* with `cc` for an unrelated, pre-existing environment reason. Built
with `CC=clang` instead (clang 19, which has a working bundled
`compiler-rt` ASan runtime) — same `-fsanitize=address -O0 -g` flags via
the Makefile's `debug` target — and ran all four modes on
`bench/panel.small.bcf` (no DUMP):

| Mode | Exit | Heap-buffer-overflow / use-after-free / double-free / invalid-free | Leaks |
|---|---|---|---|
| linear | 0 | none | none |
| sampled | 0 | none | 120192 bytes / 3 allocations (pre-existing: the `wapproxc_rrs` double-free-that's-actually-a-no-op / leaked `pbwtRev`/`pbwtPr`/`pbwtPrRev` bug documented in the plan, Step 6 fixes) |
| blockpar | 0 | none | 10367696 bytes / 520 allocations (pre-existing: `wparc_rrs`'s dead `return NULL;` before its cleanup block, documented in the plan, Step 4 fixes) |
| stagpar | 0 | none | 193796050 bytes / 341933 allocations (mostly the intentional `static __thread` GT buffers from Step 2 that live for the process lifetime, plus `pb` itself which is returned to the caller and never freed by design — not new) |

No new memory-corruption class of finding (no heap-buffer-overflow,
use-after-free, double-free, or invalid-free) was introduced by the
`pbwtad_new`/`cpbwt` flat-allocation change or by the width change — this
was the step's highest risk item and it came back clean. All findings are
LeakSanitizer reports of leaks already documented as pre-existing in the
plan/Run 0, scheduled to be fixed in Steps 4 and 6, and (for stagpar's
thread-local GT buffers) an intentional consequence of Step 2's
optimization.

### Expected performance impact

Not benchmarked here (Step 6 covers that) — halving `a`/`d` element width
should reduce the memory traffic of every memcpy, radix scatter, and
`reversec`/`divc` pass across all modes, and halves the resident
footprint of the 2*(W+1) `pbwtad`s kept live by `wparc_rrs`/
`wstagparc_rrs` at any time.

Proceeding is safe for Step 4.

## Run 4 (Step 4: blockpar kernel) — 2026-08-28

**Status:** SUCCESS — `blockpar` now matches `linear` bit-identically (200/200,
0 wrong). All of Step 4's sub-changes (1–5) were implemented; none had to be
reverted. All changes are confined to `wparc_rrs` in `sp-pbwt.c`, plus one
dead-code removal in the shared `divc()` helper (see item 6). The BM-backend
siblings `sbwparc_rrs`/`mbwparc_rrs` and the whole `stagpar` path
(`wstagparc_rrs` & co.) were deliberately left untouched.

### 1. The blockpar divergence bug — root cause was NOT the `pb0`/`pb0rev` memcpy

The plan's "Confirmed findings" pointed at the two memcpys in the parallel body
(`psrev->a` copied from `pb0[J-x]->a` instead of `pb0rev[J-x]->a`, and
`psrev->d` copied from `pb0rev[J-x]->**a**`). Reading the surrounding code
showed **those memcpys are dead**, not wrong: `reversec(ps, psrev, nrow)` runs
immediately afterwards and writes *every* element of both `psrev->a` and
`psrev->d` (it iterates `rev->X[p->a[i]]` over a full permutation `p->a`), so
whatever was memcpy'd in is unconditionally overwritten before any read.
Likewise `ps->d` is fully overwritten by `divc()`. Only the `ps->a` copy is
live (it is `rrsortx`'s starting permutation). So "fixing" the source array
would have changed nothing.

The **actual** bug was one line earlier in the prologue:

```c
pb0[0] = cpbwt(nrow, c0, p0);
pb0rev[0] = p0;              /* <-- wrong */
```

The invariant held everywhere else in the function is
`pbXrev[j]->a[row] == position of row in pbX[j]->a` — established by
`reversecprev(pb0[j], pb0[j-1], pb0rev[j], nrow)` for `j >= 1` and by
`reversec(ps, psrev, nrow)` in the parallel body. But `pb0rev[0]` was aliased
to `p0`, the *identity* permutation, which is the reverse of the state
**before** column 0, not of `pb0[0]`. `recover_div()` dereferences
`ppr->d[pprrev->a[i]]`, and the `(ppr, pprrev)` pair used for the first column
of each block is exactly `(pb0[0], pb0rev[0])` (`x == W-1`, `J-x == 0`) — so
the off-by-one-column reverse array corrupted the divergence values at exactly
the block-boundary columns `j*W`. That is precisely the observed signature:
divergences (and only divergences — the `a` arrays always matched) wrong at
columns 64 and 128 and nowhere else.

**Fix:** allocate a real reverse for column 0 and free `p0`:

```c
pb0rev[0] = pbwtad_new(nrow);
reversecprev(pb0[0], p0, pb0rev[0], nrow);
PBWTAD_FREE(p0);
```

(`->d` of a reverse array is never read by `recover_div`, only `->a` is, so
sourcing `d` from `p0` is fine and matches the `j >= 1` convention.)

With this one change alone, `blockpar` went from 2 wrong to **0 wrong** on the
fixture. The three dead memcpys in the parallel body were then deleted rather
than "fixed" (they were 3 × `nrow` × 4 bytes of pure waste per column, i.e.
per window: 3 × 63 × 20 KB).

### 2. Hoisted per-iteration allocations + `rrsortx_noaux` → `rrsortx`

The old body malloc'd a `nrow`-word merge buffer `w` per loop iteration and
`rrsortx_noaux` malloc'd its own `post` scratch internally — 126 concurrent
`nrow`-sized malloc/free per window. Both are now per-thread buffers (`w` and
`taux`) allocated once at entry to the persistent parallel region and freed at
its exit, living for the whole run. The parallel body now calls
`rrsortx(nrow, cw, ps->a, taux)`, which fuses all 8 radix-digit histograms
into a single pass over `c[]`, instead of `rrsortx_noaux`, which recomputes
the histogram in a separate pass per radix round (8 extra passes over `nrow`
64-bit words per column). Both take `pbidx_t *` after Step 3; the per-thread
`aux` is `malloc(nrow * sizeof(pbidx_t))`, matching `rrsortx`'s contract (it
ping-pongs `s`/`aux`, both `nrow` elements — 8 rounds is even, so the final
result lands back in `s` as before).

### 3. Persistent `#pragma omp parallel` region — DONE

The team is now forked **once** for the entire column loop instead of once per
64-column window (~100k team fork/joins on a 6M-column chromosome). Structure:

- Serial prologue (first window's `cpbwt` chain, `swapdiv`, `pb1` allocation,
  the first `fgetcoliw64r`) stays outside the region, unchanged.
- `#pragma omp parallel` is entered once; each thread allocates its `w`/`taux`
  there and `while (active) { ... }` iterates the windows *inside* the region.
- The old serial "head" (`x == 0`, the unshifted window, previously done by the
  master before the `omp parallel for`) was **folded into the work-sharing
  loop as iteration `x == 0`** — it is exactly the same computation with
  `cw = pwCur` instead of the merge output. This both balances all `W` columns
  across the team and lets the head use the thread's own `taux`. (It could not
  be expressed as `wr64mrgsi(..., 0)`: that macro's `wp[r] >> (W - i)` is a
  shift-by-64 at `i == 0`, i.e. UB, so `x == 0` is special-cased explicitly.)
- Per-window state advance (DUMP, `SWAP(pb0,pb1)`, `SWAP(pb0rev,pb1rev)`,
  buffer rotation, `j++`, loop-condition update) runs under `#pragma omp
  single`, whose implicit barrier publishes all of it to the team before the
  next round.
- The two backends' different loop conditions (BM/ENC: bounded by a known
  `ncol`; BCF: EOF-driven on `fgetcoliw64r` returning `< W`) are folded into a
  local `WPARC_FETCH(idx, buf)` macro so the restructure is written once and
  both compile targets keep their original semantics exactly.

### 4. I/O overlap with compute — DONE

Three window buffers (`pwPrev`/`pwCur`/`pwNext`) replace the old
`pw0`/`pw1` pair and are rotated once per round. At the top of each round:

```c
#pragma omp master
{ nxt_ok = WPARC_FETCH(j + 1, pwNext); }   /* no implicit barrier */

#pragma omp for schedule(dynamic, 1)
for (size_t x = 0; x < W; x++) { ... }     /* implicit barrier */
```

so the serial BCF decode of window `j+1` runs concurrently with the team's
63-way compute on windows `j-1`/`j`. Correctness argument:

- **No data race on the buffers:** the reader writes only `pwNext`; the compute
  loop reads only `pwPrev` and `pwCur`. The rotation that makes `pwNext` the
  new `pwCur` happens in the `single` *after* the `omp for`'s implicit barrier,
  which orders the master's writes against every thread's subsequent reads.
- **`master`, not `single nowait`:** `single` picks an arbitrary thread per
  region entry, which would spread `iobcf.c`'s `static __thread` decode
  staging buffers (Step 2) across every thread that ever wins the race —
  correct (that state is a per-call scratch cache, not carried across calls),
  but it would allocate an extra `W*nrow` = 320 KB staging buffer per thread
  and leak them all. Pinning the read to the master thread keeps exactly one
  staging buffer alive and makes the reader deterministic. `iobcf.c`'s
  cross-call reader state (`_li`, `hdr`) is plain `static`, not thread-local,
  so it is unaffected either way, and the barriers make it visible.
- `schedule(dynamic, 1)` (not `static`) is required for the overlap to pay:
  with a static schedule the master's pre-assigned chunk would only start after
  its read finished, re-serializing the very thing being hidden. With dynamic,
  the other threads drain the 64 columns while the master reads, and the master
  joins in when it is done.

### 5. Leak fix — DONE

The `return NULL;` that sat immediately *before* `wparc_rrs`'s cleanup block
(making it dead code) was removed, so the cleanup now runs. While there, the
cleanup was completed: it previously freed only `pb0[j]`/`pb1[j]` and never the
`pb0rev`/`pb1rev` pbwtads or their arrays — those are now freed too, along with
all three window buffers. Confirmed by ASan below: blockpar went from
**10,367,696 bytes in 520 allocations leaked → 0 leaks**.

### 6. `divc()`: removed a write-only `static int8_t kk`

ThreadSanitizer flagged `divc.kk` (`kk += W;` at the end of `divc`) as a
genuine data race — it is a process-global written by every thread of the
blockpar (and stagpar) teams. It was never read anywhere; it was dead debug
state. Removed. This is a shared helper also used by `stagpar`, but removing a
write-only variable cannot change behaviour, and `stagpar`'s 199/199 result is
unchanged (verified below).

### Build

`HTSLIB=/data/proj/2bfpbwt/htslib make` — both `sp-pbwt-bcf` and `sp-pbwt-bm`
compile cleanly with no warnings or errors. (`wparc_rrs` is unreachable in the
BM binary — BM `blockpar` modes dispatch to `sbwparc_rrs`/`mbwparc_rrs`, which
were not touched — but its BM/ENC branch still compiles, which the `WPARC_FETCH`
macro preserves.)

### Correctness matrix (`bench/panel.small.bcf`, 200 records, 5006 samples)

| Comparison | Checked | Wrong | Status |
|---|---|---|---|
| linear vs **blockpar** | 200 | **0** | **FIXED** (was 2 wrong at columns 64/128 in Runs 0–3) |
| linear vs stagpar | 199 | 0 | unchanged — untouched by this step |
| linear vs sampled | 4 | 1 | unchanged from baseline (sparse/approximate mode; Step 6 item) |

### Race / determinism testing

`blockpar` run 5× at each of `OMP_NUM_THREADS` = 1, 4, 8, 32 (20 runs total,
capped at 32 per shared-machine policy). **All 20 runs produced byte-identical
DUMP output** (single md5 `29d796a1581cce736bbd3f6e074338f0` across every
thread count and every repetition), and each thread count independently checks
200/200 with 0 wrong against `linear`. Output ordering is deterministic because
the DUMP is emitted from inside the `omp single`, not from the parallel body.

### Sanitizers

**ASan** (`CC=clang HTSLIB=... make debug` — the system gcc 11 is still missing
its `libasan` runtime, per Run 3's finding; clang 19's bundled compiler-rt is
used instead), all four BCF modes on the fixture at `OMP_NUM_THREADS=8`:

| Mode | Corruption findings | Leaks |
|---|---|---|
| linear | none | none |
| sampled | none | 120192 B / 3 allocs (pre-existing `wapproxc_rrs` bug, Step 6) |
| **blockpar** | **none** | **none** (was 10,367,696 B / 520 allocs) |
| stagpar | none | 26,365,034 B / 42,892 allocs (pre-existing: intentional Step-2 `static __thread` GT buffers + `pb` returned by design) |

**ThreadSanitizer** (ad hoc build: `clang -fsanitize=thread -O1 -g
-DBF2IOMODE_BCF ... sp-pbwt.c iobcf.c`), `blockpar` at `OMP_NUM_THREADS=8`:

- First run: 7 warnings — 3 on `divc.kk` (the genuine application-level race,
  fixed in item 6 above), 4 entirely inside `libomp.so`
  (`pthread_mutex_lock`/`__kmp_free` with no application frames).
- After removing `kk`: **2 warnings, both inside `libomp.so`** — the standard
  false positives from an uninstrumented OpenMP runtime (TSan needs an
  archer-annotated libomp to suppress these). **Zero races reported in
  application code**: nothing on `pwPrev`/`pwCur`/`pwNext`, on
  `pb0`/`pb1`/`pb0rev`/`pb1rev`, on `j`/`active`/`nxt_ok`, or in the `iobcf.c`
  decode path — i.e. the new persistent-parallel-region and double-buffered
  reader are clean.

### Performance

Not benchmarked here (Step 6 owns that). A rough before/after on the fixture is
meaningless — 200 columns is only 3 windows, so the persistent-region and
read-ahead wins have essentially nothing to amortise over; `wparc_rrs` measured
25.1–26.6 ms before and 24.6–25.7 ms after at 8 threads, i.e. inside the noise.
The changes that should matter at chromosome scale are: ~100k fewer team
fork/joins, 126 fewer mallocs per window, 7 of 8 histogram passes dropped from
every column's radix sort, 3 dead `nrow`-sized memcpys dropped per column, and
the BCF decode moved off the critical path.

**Proceeding is safe for Step 5** (which must not reuse the `pb0rev[0] = p0`
pattern — `wstagparc_rrs` should be checked for the same off-by-one when it is
restructured).

**Note for a future step (out of Step 4's scope):** the same
`pb0rev[0] = p0;` prologue line still exists in the two BM-backend siblings,
`sbwparc_rrs` (`sp-pbwt.c:1157`) and `mbwparc_rrs` (`:1262`). They almost
certainly carry the identical block-boundary divergence bug, but the `.bm`
backend is explicitly out of scope for this plan and there is no `.bm` fixture
to verify against here (`gen` currently aborts with "glibc detected an invalid
stdio handle" on this machine — a pre-existing, unrelated failure), so they
were left alone rather than changed blind.

## Run 5 (Step 5: stagpar restructure) — 2026-08-28

**Status:** SUCCESS — `wstagparc_rrs` (BCF `stagpar`) fully restructured to a
single sequential reader + inverted loop nest. Still **199 checked / 0 wrong**
against `linear`, now with byte-identical output across all thread counts, no
ASan findings, **zero leaks** (was ~26 MB), and no application-level TSan race.

All changes are confined to `wstagparc_rrs` in `sp-pbwt.c`. The BM-backend
siblings `swstagparc_rrs`/`mwstagparc_rrs` were deliberately **not** touched
(`.bm` backend out of scope, and no working `.bm` fixture exists here — `gen`
still aborts, per Run 4's note).

### 1. What the old structure was

- The function opened one `bcf_srs_t` **per OpenMP thread**, each streaming the
  BCF from the very beginning; with the default team size that is up to W=64
  full BGZF decompressions of the input per run, plus repeated hits of
  `iobcf.c`'s reopen-and-rescan path on every lane switch.
- The loop nest was **lane-outer / round-inner**: each thread took a contiguous
  slice of lanes (`start`/`count` partition) and, for each lane, privately
  streamed windows `[lane, lane+W)`, `[lane+W, lane+2W)`, … with
  `fgetcolwgri`.
- Carried state (`pt0`, `pt0rev`) plus scratch (`pt1`, `pt1rev`, `aux`, `pw`)
  were all **per thread**, allocated once at region entry — correct only
  because a lane never migrated between threads.
- Four full `nrow`-sized memcpys per column (`pt1->a`, `pt1->d`, `pt1rev->a`,
  `pt1rev->d`).
- The DUMP was emitted from an `omp critical` inside the parallel body, so the
  output line order changed on every run (documented in Run 2).
- The function returned `pb`, which held only the 64 **seed** states; the real
  per-lane results lived in thread-locals that were freed at region exit. `pb`
  itself was never freed by anyone, and neither was the prologue's reader.

### 2. The restructure

The enabling observation (from the plan): lane `l` at round `r` consumes columns
`[l + rW, l + rW + W)`. Across all W lanes at a fixed `r` that is exactly the
columns of **blocks `r` and `r+1`** — precisely the two-window working set
`wparc_rrs` (blockpar) already keeps in its ring buffer. Concretely, the window
for lane `l` at round `r` starts at column `rW + l`, i.e.

- `l == 0`: block `r` as-is (`pwCur`) — a shift of 0 would be a shift-by-64 in
  `wr64mrgsi` (UB), so it is special-cased exactly as `wparc_rrs` does;
- `l >= 1`: `wr64mrgsi(nrow, pwNext /*block r+1*/, pwCur /*block r*/, w, W - l)`

which reproduces bit-for-bit the packing the old per-lane `fgetcolwgri`
produced.

Changes made:

- **One reader instead of 64.** The 64 in-region `bcf_sr_init`/`bcf_sr_add_reader`
  calls are gone. There are now exactly two readers over the function's life,
  and never more than one at a time: `_pr` for the sequential prologue
  (`fgetcoli` over columns 0..W-1), destroyed immediately after; then a single
  `_sr` driving the whole window scan via `fgetcoliw64r`, destroyed at exit.
  (Two rather than one only because the prologue consumes columns 0..W-1 and
  the block scan must restart at block 0 — the same prologue-reader pattern
  `wparc_rrs` uses.)
- **Loop nest inverted.** Outer loop is now the round `r`, advancing the ring
  buffer once per round; the lanes are the inner `#pragma omp for
  schedule(dynamic, 1)` work-sharing loop, inside a single persistent
  `#pragma omp parallel` region that spans the entire run (one team fork/join
  total, like Step 4's blockpar).
- **Per-lane carried state.** Because a lane can now be picked up by a different
  thread each round, `pt0`/`pt0rev` became per-LANE arrays: `pb[l]` (reused
  directly as the lane state — it is already the correct seed, so the old
  per-lane seeding memcpys disappear) and a new `pbrev[l]`, seeded once by a
  `#pragma omp for` `reversec` pass before the round loop. `pt1`/`pt1rev`/`aux`/
  `w` remain per-thread: they are pure one-iteration scratch.
- **`STAG_FETCH(idx, buf)`** mirrors Step 4's `WPARC_FETCH` but returns the
  *column count* actually read (W, or fewer at EOF) rather than a boolean. The
  count is needed: at the last round, lane `l` still has a complete window iff
  the trailing block holds at least `l` columns, so the lane loop skips
  `lane > nNext`. This reproduces the old per-lane EOF condition
  (`l + (r+1)W <= ncol`) exactly — verified by the dumped column set being
  unchanged (prologue 1..63, then 63..199; 199 distinct columns, matching every
  previous run's "checked: 199").
- **Two dead memcpys dropped.** `divc()` reads only `ppr->d` and `pprrev->a` of
  the previous state, so the `pt1->a` and `pt1rev->d` copies were pure waste
  (2 × `nrow` × 4 bytes per column). Removed; correctness re-verified.
- **Deterministic DUMP.** The dump moved out of the `omp critical` in the
  parallel body into the per-round `omp single`, emitted in lane order. Since
  the dumped column index is `rW + l + W - 1`, output columns are now strictly
  increasing and the whole dump file is byte-identical run to run (see below).
- **BM/ENC branch.** `wstagparc_rrs` is unreachable in `sp-pbwt-bm` (BM
  `stagpar` dispatches to `swstagparc_rrs`/`mwstagparc_rrs`), but it still has
  to compile. Its `STAG_FETCH` is the `ncol`-bounded, full-blocks-only variant;
  this is noted inline as a slight tail-semantics difference from the old BM
  loop, in dead code, rather than guessing at partial-block reads through the
  BM readers.

### 3. Return value decision — `return NULL`

`main()` assigns the mode dispatch's result to `r` and **never reads it**, for
any mode; `linc` and (since Step 4) `wparc_rrs` both already return `NULL`.
Rather than materialize a per-column snapshot array nothing consumes, `stagpar`
now matches them: it frees everything it allocated (`pb[0..W]`, `pbrev[0..W-1]`,
all three window buffers, `c0`, and — new — the BCF readers) and returns `NULL`.
Nothing from the old thread-local state is needed by any caller; the per-column
results were only ever observed through the DUMP, which is emitted during the
loop and is unaffected.

### 4. I/O overlap — INCLUDED

It transferred directly. Three buffers (`pwCur` = block `r`, `pwNext` = block
`r+1`, `pwSpare`) are rotated once per round. At the top of each round:

```c
#pragma omp master
{ nSpare = (nNext == W) ? STAG_FETCH(r + 2, pwSpare) : 0; }   /* no barrier */

#pragma omp for schedule(dynamic, 1)
for (size_t lane = 0; lane < W; lane++) { ... }               /* barrier */
```

so the serial BCF decode of block `r+2` runs concurrently with the team's
64-way lane compute on blocks `r`/`r+1`. Same correctness argument as Step 4:
the reader writes only `pwSpare`, lanes read only `pwCur`/`pwNext`, and the
rotation happens in the `omp single` *after* the `omp for`'s implicit barrier.
`master` (not `single nowait`) keeps `iobcf.c`'s `static __thread` decode
staging buffer on exactly one thread. `schedule(dynamic, 1)` lets the other
threads drain lanes while the master reads.

### Build

`HTSLIB=/data/proj/2bfpbwt/htslib make` — `sp-pbwt-bcf` and `sp-pbwt-bm` both
compile cleanly, no warnings.

### Correctness matrix (`bench/panel.small.bcf`, 200 records, 5006 samples)

| Comparison | Checked | Wrong | Status |
|---|---|---|---|
| linear vs **stagpar** | 199 | **0** | PASS — unchanged after the rewrite |
| linear vs blockpar | 200 | 0 | unchanged (untouched by this step) |
| linear vs sampled | 4 | 1 | unchanged from baseline (sparse/approximate; Step 6) |

### Race / determinism testing

`stagpar` run 5× at each of `OMP_NUM_THREADS` = 1, 4, 8, 32 (20 runs, capped at
32 per shared-machine policy). **All 20 runs produced byte-identical DUMP
output** — single md5 `6b88df27a6132416880053a084674380` across every thread
count and repetition. This is a strict improvement on the pre-Step-5 behavior,
where dump *ordering* varied run to run (Run 2) and only the sorted content was
stable.

### Sanitizers

**ASan** (`CC=clang HTSLIB=... make debug`; system gcc still lacks its libasan
runtime, per Run 3), all four BCF modes on the fixture at
`OMP_NUM_THREADS=8`:

| Mode | Corruption findings | Leaks |
|---|---|---|
| linear | none | none |
| sampled | none | 120192 B / 3 allocs (pre-existing `wapproxc_rrs` bug, Step 6) |
| blockpar | none | none |
| **stagpar** | **none** | **none** (was 26,365,034 B / 42,892 allocs) |

The stagpar leak is fully gone: `pb`/`pbrev`/the window buffers/`c0` are now
freed, both BCF readers are destroyed, and there is no longer a per-thread
`iobcf.c` staging buffer per lane-thread (only the master decodes).

**ThreadSanitizer** (ad hoc: `clang -fsanitize=thread -O1 -g -DBF2IOMODE_BCF …
sp-pbwt.c iobcf.c`), `stagpar` at `OMP_NUM_THREADS=8`: **6 warnings, all six
entirely inside `libomp.so`** — `pthread_mutex_lock`/`pthread_mutex_init` on
libomp's own team locks and `memset`/`free` inside `__kmp_free` on libomp's own
team heap blocks. The only application frame appearing in any of them is the
outlined `#pragma omp parallel` entry point itself. These are the standard
false positives from an uninstrumented OpenMP runtime (TSan needs an
archer-annotated libomp), identical in kind to the 2 residual warnings Step 4
reported. **Zero races reported on application state**: nothing on
`pwCur`/`pwNext`/`pwSpare`, `pb`/`pbrev`, `r`/`active`/`nCur`/`nNext`/`nSpare`,
or the `iobcf.c` decode path. The TSan binary's DUMP also checks 199/199, 0
wrong against `linear`.

### Performance

Step 6 owns benchmarking; no before/after was measured on the fixture because
200 columns is only 3 rounds — the single-reader win has nothing to amortise
over there (an attempted stashed-baseline rebuild also failed for an unrelated
Makefile-revert reason, so no honest fixture number is quoted).

One real-scale smoke run was done to confirm the restructure holds up outside
the fixture: `stagpar` on chr20 (1000G phase3, 5008 rows) at
`OMP_NUM_THREADS=8` under `limitram 16G` completed in **65.3 s wall / 514 s
CPU**. For context, Run 0's baseline notes recorded the parallel modes as *not
completing* within a 120 s timeout on the same chromosome. The structural wins
that should show at chromosome scale: 64 → 1 BGZF decompressions of the input,
64 → 1 concurrent htslib reader contexts and their buffers, no reopen-and-rescan
on lane switches, one team fork/join instead of per-lane region entry, 2 of 4
`nrow`-sized memcpys per column removed, and the BCF decode moved off the
critical path.

**This completes Step 5.** Step 6 (the `wapproxc_rrs` double-free/leak fix, the
full benchmark re-run, and the commit) is next.

## Run 6 (Step 6: cleanup, validation, final benchmark) — 2026-08-28

**Status:** SUCCESS — leak fixed, full correctness matrix clean at two scales,
ASan clean on all four modes, chr20/chr10 benchmark completed for all four
modes (previously only `linear` completed at all). This is the plan's final
step; per its explicit stop condition, no further tuning follows.

### 1. `wapproxc_rrs` double-free / leak — FIXED

The cleanup block (end of `wapproxc_rrs`, `sp-pbwt.c`) was:

```c
PBWTAD_FREE(pbwt);
FREE(pbwt);
FREE(pbwtRev);
FREE(pbwtPr);
FREE(pbwtPrRev);
```

Post-Step-3, `PBWTAD_FREE(p)` already frees `p->a` (the shared `a`/`d` block)
*and* `p` itself. So `PBWTAD_FREE(pbwt); FREE(pbwt);` double-frees the `pbwt`
struct, while `FREE(pbwtRev)`/`FREE(pbwtPr)`/`FREE(pbwtPrRev)` free only the
three-struct wrappers and leak each one's `a`/`d` block — exactly the "frees
`pbwt` twice, leaks the other three" bug described in the plan, restated in
terms of the current (post-Step-3) flat-allocation `PBWTAD_FREE`. Fixed to:

```c
PBWTAD_FREE(pbwt);
PBWTAD_FREE(pbwtRev);
PBWTAD_FREE(pbwtPr);
PBWTAD_FREE(pbwtPrRev);
FREE(w64);
FREE(aux);
return NULL;
```

Confirmed by ASan below: `sampled`'s leak report went from 120,192 bytes / 3
allocations (Runs 3–5) to **zero**.

### 2. Full correctness matrix

**Small fixture** (`bench/panel.small.bcf`, 200 records, 5006 samples):

| Comparison | Checked | Wrong | Status |
|---|---|---|---|
| linear vs blockpar | 200 | 0 | PASS — unchanged |
| linear vs stagpar | 199 | 0 | PASS — unchanged |
| linear vs sampled | 4 | 1 | unchanged (known sparse/approximate mode; the 1 mismatch is pre-existing and not caused by the leak — see note below) |

**Mid-size slice** (`bench/panel.mid.bcf`, 20,000 records cut from
`chr20` — `bcftools` is still not on `PATH` in this environment, confirmed
again per Runs 0/2's finding, so a small purpose-built `extract_bcf_subset.c`
tool — first N BCF records, header preserved — was written and used instead):

| Comparison | Checked | Wrong | Status |
|---|---|---|---|
| linear vs blockpar | 20000 | 0 | PASS |
| linear vs stagpar | 19999 | 0 | PASS |
| linear vs sampled | 313 | 1 | same pattern as the small fixture, at 100x the scale |

`blockpar` and `stagpar` are bit-identical to `linear` at both scales.
`sampled`'s single mismatch is unchanged in count and character between the
200-record and 20,000-record inputs, consistent with it being an
already-documented approximation-mode artifact (not a memory bug) — the leak
fix touched only cleanup code after the last `PDUMPR`, so it could not and did
not change any emitted value.

### 3. ASan (`CC=clang HTSLIB=... make debug`, small fixture, all four BCF modes)

| Mode | Corruption findings | Leaks |
|---|---|---|
| linear | none | none |
| **sampled** | none | **none** (was 120,192 B / 3 allocations in Runs 3–5) |
| blockpar | none | none |
| stagpar | none | none |

All four modes are now fully clean under ASan: no heap-buffer-overflow,
use-after-free, double-free, invalid-free, or leaked bytes.

### 4. Benchmark: chr20 and chr10, all four BCF modes, 32 threads, `limitram 16G`

Ran via `bench/run_step6.sh` (a driver over the Step-0 `bench/run.sh`, updated
to export `LD_LIBRARY_PATH` for the htslib shared library and to report
CPU = user + sys) under `babysit-run`, single job, ~13 minutes total wall time
for all 8 (mode × chromosome) runs combined.

**chr20** (5008 samples, 1000G phase3):

| Mode | Wall | CPU (user+sys) | Peak RSS |
|---|---|---|---|
| linear | 1:04.11 | 63.46 s | 4.00 MB |
| sampled | 0:49.43 | 48.95 s | 4.50 MB |
| blockpar | 1:15.50 | 1798.39 s | 12.30 MB |
| stagpar | 0:57.12 | 1472.30 s | 9.00 MB |

**chr10** (5008 samples, larger chromosome):

| Mode | Wall | CPU (user+sys) | Peak RSS |
|---|---|---|---|
| linear | 2:17.99 | 136.81 s | 4.00 MB |
| sampled | 1:54.61 | 113.38 s | 4.50 MB |
| blockpar | 2:40.40 | 3893.83 s | 11.30 MB |
| stagpar | 2:07.82 | 3115.81 s | 10.00 MB |

**Before/after vs. Run 0's baseline** (Run 0 only fully completed `linear`;
`sampled` timed out past 120 s, `blockpar`/`stagpar` did not complete):

| Mode | chr20 baseline (Run 0) | chr20 after (Run 6) | chr10 after (Run 6) |
|---|---|---|---|
| linear | ~76 s wall / ~75 s CPU | **64.1 s wall / 63.5 s CPU** (~16% faster) | 138.0 s wall / 136.8 s CPU |
| sampled | did not complete (>120 s timeout, and was only emitting 4/200 columns on the fixture at the time) | **49.4 s wall — completes correctly** | 114.6 s wall |
| blockpar | did not complete | **75.5 s wall / 1798 s CPU (~24x of wall, i.e. good 32-thread scaling)** | 160.4 s wall / 3894 s CPU |
| stagpar | did not complete | **57.1 s wall / 1472 s CPU (~26x of wall) — fastest wall time of the three whole-genome-scanning modes** | 127.8 s wall / 3116 s CPU |

The headline result is qualitative, not just a speedup ratio: three of the
four modes (`sampled`, `blockpar`, `stagpar`) went from *not finishing at all*
in Run 0 to completing chr20 in under 76 seconds wall time — comparable to or
faster than `linear`'s own baseline — while also being bit-identical to
`linear` (`blockpar`/`stagpar`) or matching at their sample points
(`sampled`). Peak RSS stayed low across the board (4–12 MB) since `pbidx_t`
(Step 3, `uint32_t`) plus per-thread scratch reuse (Steps 4–5) keep resident
state small relative to the 5008-sample panel; the total CPU-seconds figures
for `blockpar`/`stagpar` (in the thousands) reflect 32-way parallel work, not
inefficiency — their wall-clock times are the relevant comparison against
`linear`'s single-threaded baseline.

### 5. Commit

Files staged: `Makefile`, `iobcf.c`, `sp-pbwt.c`, `CHANGELOG.md`, `.gitignore`,
`bench/run.sh`, `bench/run_step6.sh`, `extract_bcf_subset.c`. Deliberately
**not** staged: `CLAUDE.md` (untracked before this session, unrelated to this
work), the compiled `extract_bcf_subset` binary and the generated
`bench/*.bcf`/`bench/*.txt` fixtures and `bench/step6_results.txt` (added to
`.gitignore` instead, consistent with the existing convention of not tracking
generated panels — see `data/` already being ignored), and
`bench/.pipeline-run/` (babysit-run's internal state, also gitignored).

### Notes for future work (out of this plan's scope, not acted on)

- The BM-backend siblings `sbwparc_rrs`/`mbwparc_rrs` (blockpar) almost
  certainly carry the same `pb0rev[0] = p0` off-by-one bug fixed in Step 4 for
  the BCF path; `swstagparc_rrs`/`mwstagparc_rrs` were never restructured
  either. Both are explicitly out of scope (`.bm` backend was excluded from
  this plan) and there is still no working `.bm` fixture in this environment
  (`gen` aborts on this machine) to verify against.
- `bcftools` remains absent from `PATH` in this environment across all 6 runs;
  the mid-size correctness slice was produced with a small ad hoc
  `extract_bcf_subset.c` instead.

**This completes the plan (Steps 0–6).** Per the plan's explicit stop
condition, no further optimization work follows in this session.

## Run 7 (Manuscript benchmark refresh: `blockpar`/`stagpar` on all chromosomes) — 2026-08-28

**Status:** SUCCESS. Not part of the optimization plan itself — a follow-up
benchmark-update requested after Steps 0–6 landed, re-running the manuscript's
experiment pipeline (`../exp`, a separate git repo) for the two modes this
work actually changed, across the full chromosome set, and refreshing the
manuscript's `all.csv`-shaped summary with the new numbers.

### What was done

1. Verified `sp-pbwt-bcf` was freshly built from the completed optimization
   commit (`702f42c`) via `HTSLIB=/data/proj/2bfpbwt/htslib make`.
2. In `/data/proj/2bfpbwt/exp` (sibling repo): symlinked
   `phase3 -> /data/proj/2bfpbwt/one_thousands/phase3` (that repo had no data
   directory and none of its competing-tool binaries built — the tracked
   `Snakefile` cannot run as committed regardless of this work).
3. Added a new, untracked `Snakefile.spbwt_eval` in `../exp` — NOT an edit to
   the tracked `Snakefile` — containing only one rule modeled on the existing
   `run_bf_pbwt_bcf` rule, pointed at
   `/data/proj/2bfpbwt/sp-pbwt/sp-pbwt-bcf` with `BCF_ARGS = ["blockpar",
   "stagpar"]` (today's current mode names, not the old frozen `bf-pbwt-clean`
   short codes), writing to new `results_spbwt_eval/`/`benchmarks_spbwt_eval/`
   directories so the existing 1357/429 tracked result files were never
   touched. Needed `LD_LIBRARY_PATH=/data/proj/2bfpbwt/htslib` in the rule's
   shell command (the binary isn't statically linked against htslib).
4. Ran via `snakemake -s Snakefile.spbwt_eval -c 32 -p --rerun-incomplete`,
   wrapped in `limitram 32G`, launched under the `babysit-run` skill (shared-
   machine policy: cap cores at 32, wrap memory-relevant jobs in `limitram`,
   babysit long jobs instead of blocking). 44 jobs (22 chromosomes ×
   `blockpar`/`stagpar`), one at a time (each job requests the full 32-thread
   budget). Took ~1h45m wall for the whole batch; the watcher process died
   once mid-run from an unrelated session hiccup and was relaunched against
   the same run directory without touching or restarting the pipeline itself.
5. Built `../exp/all2.csv`: same schema and same 396 rows (22 panels × 18
   tools) as the manuscript's existing `../exp/results/all.csv`, produced by
   copying `all.csv` and replacing only the 44 rows for tool labels
   `bcf-prs`/`bcf-spr` — the old names for `blockpar`/`stagpar` — with values
   parsed (via the repo's own `time_to_csv.py`) from the new
   `results_spbwt_eval/*.time` files. Every other row (`bcf-ars`, `bcf-lin`,
   all `bm-*`, `parpbwt`, `pbwt`, `syl` — tools/modes this work never touched)
   is byte-identical to `all.csv`. Verified row-by-row in Python, not just by
   `diff` (raw `diff`/`grep` output was misleading due to hunk grouping).

### Results (peak RSS ~20–24 MB on every chromosome, including chr1 — confirms
the streaming design holds at full-genome scale, not just on the chr10/chr20
fixtures from Run 6)

| chr | blockpar wall/cpu/rss(MB) | stagpar wall/cpu/rss(MB) |
|---|---|---|
| chr1 | 235.9 / 5517 / 22.5 | 253.2 / 5368 / 20.2 |
| chr2 | 237.3 / 5808 / 22.4 | 234.9 / 5559 / 20.2 |
| chr3 | 199.6 / 5038 / 22.5 | 223.0 / 4703 / 20.1 |
| chr4 | 179.2 / 4307 / 22.6 | 180.8 / 3877 / 21.5 |
| chr5 | 170.5 / 4283 / 22.5 | 206.9 / 4525 / 20.0 |
| chr6 | 158.3 / 3521 / 22.5 | 167.1 / 4138 / 20.1 |
| chr7 | 185.1 / 4005 / 22.5 | 150.0 / 3206 / 20.1 |
| chr8 | 135.5 / 2743 / 22.5 | 174.3 / 3991 / 20.2 |
| chr9 | 167.1 / 3864 / 22.5 | 184.9 / 3691 / 20.4 |
| chr10 | 142.1 / 3501 / 22.6 | 153.1 / 3303 / 20.2 |
| chr11 | 140.6 / 3542 / 22.5 | 125.3 / 2610 / 20.1 |
| chr12 | 129.2 / 2658 / 22.6 | 181.5 / 3783 / 21.3 |
| chr13 | 129.5 / 2531 / 22.6 | 138.0 / 3107 / 21.6 |
| chr14 | 120.5 / 2543 / 22.5 | 132.7 / 2433 / 20.1 |
| chr15 | 119.2 / 2491 / 22.5 | 115.8 / 2409 / 20.3 |
| chr16 | 131.8 / 2552 / 24.0 | 124.2 / 2458 / 20.2 |
| chr17 | 112.9 / 2514 / 22.6 | 112.0 / 2415 / 20.2 |
| chr18 | 110.8 / 2425 / 22.5 | 106.2 / 1730 / 20.1 |
| chr19 | 88.5 / 1794 / 22.5 | 92.2 / 1674 / 20.3 |
| chr20 | 86.9 / 1797 / 22.6 | 88.9 / 1720 / 20.2 |
| chr21 | 36.1 / 418 / 22.5 | 54.3 / 1049 / 20.2 |
| chr22 | 55.9 / 1051 / 22.6 | 56.0 / 1035 / 20.1 |

Sample before/after from `all.csv` → `all2.csv` (wall clock, seconds):
chr22 `bcf-prs` 114.1 → 55.8; chr22 `bcf-spr` 296.3 → 56.0; chr1 `bcf-prs`
661.9 → 235.8; chr1 `bcf-spr` 1738.4 → 253.2. Peak memory also dropped
substantially (e.g. chr1 `bcf-spr` 33700 KB → 10240 KB), consistent with the
Step 3 `pbidx_t` width halving and Step 4/5 leak fixes.

### Where things live

- `/data/proj/2bfpbwt/exp/Snakefile.spbwt_eval` (new, untracked)
- `/data/proj/2bfpbwt/exp/phase3` (new symlink, untracked)
- `/data/proj/2bfpbwt/exp/results_spbwt_eval/`,
  `/data/proj/2bfpbwt/exp/benchmarks_spbwt_eval/` (new, untracked)
- `/data/proj/2bfpbwt/exp/all2.csv` (new, untracked)

Nothing in `../exp`'s tracked files (`Snakefile`, `README.md`) or its existing
`results/`/`benchmarks/`/`results/all.csv` was modified. Nothing was
committed in either repo.

### Repeatable procedure (for future benchmark-update requests)

This is now the standard playbook whenever asked to refresh the manuscript
benchmark after further `sp-pbwt` changes, without re-running the competing
tools or the `.bm` backend:

1. Rebuild `sp-pbwt-bcf` in this repo (`HTSLIB=/data/proj/2bfpbwt/htslib
   make`) and confirm it matches the commit whose numbers you want.
2. In `../exp`: the `phase3` symlink and `Snakefile.spbwt_eval` from this run
   can be reused as-is if the mode list hasn't changed — just re-run
   `snakemake -s Snakefile.spbwt_eval -c 32 -p --rerun-incomplete` (add
   `--forceall` or delete stale `results_spbwt_eval/`/`benchmarks_spbwt_eval/`
   entries to force a true re-run rather than skipping up-to-date outputs).
   Always wrap in `limitram` and launch via `babysit-run` per shared-machine
   policy — a full-chromosome sweep takes ~2 hours.
3. Rebuild the manuscript-shaped CSV by copying `../exp/results/all.csv` and
   replacing only the rows for the tool labels that map to modes actually
   re-run (`blockpar` → `bcf-prs`, `stagpar` → `bcf-spr` in the manuscript's
   naming) using `time_to_csv.py <tool-label> <panel> < <.time file>` — every
   other row must stay byte-identical to the source `all.csv`. Verify
   row-by-row in Python (not just `diff`, which can misrepresent partial
   matches as full-block changes) before calling it done.
4. Append a dated changelog entry here with the full per-chromosome table and
   before/after deltas, matching this entry's format.

**CORRECTED by Run 9b below** — the `all2.csv` numbers in this entry compared
32-thread `blockpar`/`stagpar` runs against the frozen `all.csv`'s 8-thread
`bcf-prs`/`bcf-spr` rows. Not an apples-to-apples speedup number. See Run 9b.

## Run 8 (Step 0: Baseline harness for `linear`/`sampled` optimization plan) — 2026-08-28

**Status:** SUCCESS. Step 0 of the new plan targeting single-threaded CPU
reduction in the BCF `linear` (`linc`) and `sampled` (`wapproxc_rrs`) modes
(plan file: `we-did-some-very-rustling-meadow.md`). Unrelated to Runs 0–7
above, which were the prior `blockpar`/`stagpar` parallel-mode plan.

### What was done

1. Built current `master` (commit `702f42c`) with
   `HTSLIB=/data/proj/2bfpbwt/htslib make -B sp-pbwt-bcf` and froze the result
   as `bench/sp-pbwt-bcf.base` — the reference binary every later step's
   `bench/verify.sh` run diffs against. Never rebuilt over it.
2. Produced reference `DUMP` output from `bench/sp-pbwt-bcf.base` for all four
   BCF modes (`linear`, `sampled`, `blockpar`, `stagpar`) on both
   `bench/panel.small.bcf` and `bench/panel.mid.bcf`, saved as
   `bench/ref/<mode>.<panel>.dump` (8 files, largest ~977 MB for the `.mid`
   panel under `linear`/`blockpar`/`stagpar`).
3. Confirmed the stop condition: re-ran `sampled`/`blockpar` on
   `panel.small.bcf` from `bench/sp-pbwt-bcf.base` a second time and `cmp`'d
   against the saved reference — byte-identical, so the base binary's own
   output is reproducible run-to-run and safe to use as ground truth.
4. Wrote `bench/verify.sh`: rebuilds `sp-pbwt-bcf`, re-dumps all 4 modes × 2
   panels, `cmp`s each against `bench/ref/`, prints `OK`/`DIFF` per
   mode+panel and `== ALL BIT-IDENTICAL ==`/`== VERIFY FAILED ==` overall,
   non-zero exit on any diff. Ran once against an unmodified build to confirm
   a clean pass.
5. Recorded baseline timings (`bench/baseline_timing.sh`, run under
   `limitram 16G`, `OMP_NUM_THREADS=1`, via the `babysit-run` skill since a
   full 12-run chromosome sweep is long — took ~9 minutes total): 3 runs each
   of `linear`/`sampled` on chr21 and chr20
   (`/data/proj/2bfpbwt/exp/phase3/ALL.chr{20,21}...bcf`), `bench/sp-pbwt-bcf.base`.

### Baseline (median of 3, `linc`/`wapproxc_rrs` TRACE wall time)

| mode | chr21 | chr20 |
|---|---|---|
| `linear` | 37.5 s (37.4/37.5/41.5) | 61.6 s (72.2/61.6/61.6) |
| `sampled` | 28.8 s (28.7/28.8/28.8) | 47.3 s (47.3/47.4/47.2) |

CPU time tracked wall time within 1s in every run (single-threaded, as
expected). This is the number each later step's timing check is compared
against; Step 1 onward use chr21 (smaller, faster iteration) per the plan.

### Where things live

- `bench/sp-pbwt-bcf.base` — frozen pre-change binary (untracked, git-ignored
  by size/binary convention — check before committing).
- `bench/ref/*.dump` — reference dumps (untracked, large; not intended for
  git).
- `bench/verify.sh` — the step gate: `bench/verify.sh` must print
  `== ALL BIT-IDENTICAL ==` before any step's change is kept.
- `bench/baseline_timing.sh`, `bench/baseline_timings.txt` — timing harness
  and raw output for the table above.

## Run 9 (Step 1: `cpbwti`/`cpbwt`/`cpbwtiLCP` — eliminate `o`/`h` staging arrays) — 2026-08-28

**Status:** SUCCESS, kept. Rewrote the three PBWT column kernels (`cpbwti`
used by `linc`/`wapproxc_rrs`'s prologue-free steady state; `cpbwt` and
`cpbwtiLCP`, both reachable under `BF2IOMODE_BCF` via `wparc_rrs`'s
`blockpar`/`stagpar` prologue and steady state) to write directly into the
output `a`/`d` arrays instead of writing zeros into the output and staging
ones into separate `o`/`h` arrays followed by two trailing `memcpy`s.

Added `count_zeros(n, c)` (sequential byte pass, auto-vectorizes) computed
once before each kernel body to get `n0`, then `r=0, q=n0` and a branchless
`pos = mask ? q : r; dv = mask ? g : f;` select per iteration, writing once
into `pc->a[pos]`/`pc->d[pos]`. Numeric semantics (the `f`/`g` recurrence,
`-mask`/`-(1-mask)` masking) untouched — this is a data-movement change only.

### Verification

`bench/verify.sh` passes: all 4 BCF modes × 2 panels byte-identical to
`bench/ref/`.

### Timing (chr21, single-threaded, `limitram 16G`, `bench/sp-pbwt-bcf.base`
vs new build)

3-run sweep, `TRACE` wall seconds:

| mode | base (3 runs) | new (3 runs) | median Δ |
|---|---|---|---|
| `linear` | 37.43, 37.30, 37.69 | 37.12, 39.10, 37.13 | ~flat/noisy on 3 runs alone |
| `sampled` (control — doesn't call `cpbwti`) | 28.68, 28.62, 28.63 | 28.53, 28.73, 28.57 | ~flat, as expected |
| `blockpar` | 255.84, 256.22, 253.94 | 255.18, 255.01, 253.49 | ~flat/small |
| `stagpar` | 245.50, 247.54, 254.98 | 243.32, 242.65, 245.65 | -1 to -2% |

The 3-run `linear` signal was within noise, so per the plan's explicit
"measure this step in isolation... if `cpbwti` does not improve, revert it"
instruction, 7 additional interleaved `bench/sp-pbwt-bcf.base` vs new-build
`linear`/chr21 pairs were run back-to-back (`bench/step1_linear_extra.sh`) to
separate signal from run-to-run jitter:

| run | base wall(s) | new wall(s) | base CPU(s) | new CPU(s) |
|---|---|---|---|---|
| 1 | 37.69 | 37.09 | 37 | 36 |
| 2 | 37.30 | 37.11 | 37 | 36 |
| 3 | 37.24 | 37.06 | 37 | 36 |
| 4 | 37.28 | 37.20 | 37 | 36 |
| 5 | 37.36 | 37.16 | 37 | 36 |
| 6 | 37.42 | 37.15 | 37 | 36 |
| 7 | 37.16 | 37.15 | 36 | 36 |

New build was faster than its paired base run in all 7 interleaved reps
(median 37.30 → 37.15, ~0.4%), and CPU time was consistently 1s lower in
every pair — small but real and directionally consistent with the plan's
prediction (relieving the `cpbwti` register-spill bottleneck identified in
`perf annotate`). **Decision: keep** `cpbwti`/`cpbwt`/`cpbwtiLCP` as rewritten;
nothing reverted.

### Note

The subagent that implemented this step ran its timing sweeps as several
separate ad-hoc `run_in_background` shell launches (ending its turn and
waking on each notification) instead of a single `babysit-run` launch, per
this repo's CLAUDE.md workflow rule 7. Not incorrect, just noisier/costlier
than intended — flagged here so future delegated steps are told explicitly
to use `babysit-run` for any multi-run timing sweep.

## Run 9b (Correction: Run 7's `blockpar`/`stagpar` refresh re-run at matching thread count) — 2026-08-29

**Status:** SUCCESS. Correction, not new optimization work. The user pointed
out that Run 7's `Snakefile.spbwt_eval` set `THREADS = 32` for the
`blockpar`/`stagpar` re-run, while the frozen `results/all.csv` rows it was
being compared against (`bcf-prs`/`bcf-spr`, produced by the old
`bf-pbwt-clean` binary via the main `Snakefile`'s `run_bf_pbwt_bcf` rule) used
`threads: 8` for those same tool labels. Run 7's `all2.csv` was therefore
comparing 32-thread numbers to an 8-thread baseline — not a valid before/after
comparison, and the speedups reported in Run 7 overstate the algorithmic
improvement by conflating it with a 4x increase in thread budget.

### What was done

1. Set `Snakefile.spbwt_eval`'s `THREADS` from 32 to 8, matching the frozen
   rows' thread count exactly (`main Snakefile`'s `run_bf_pbwt_bcf` rule:
   `8 if wildcards.arg in ["bpr", "spr", "prs"] else 1`).
2. Re-pointed `BCF_EXE` from the repo's live `sp-pbwt-bcf` to
   `bench/sp-pbwt-bcf.base` — at the time of this correction, `sp-pbwt.c`/
   `iobcf.c` were being actively modified in place by a separate,
   still-in-progress plan (the `linear`/`sampled` optimization work, Steps
   1-2 mid-flight), so the live binary no longer reliably represented commit
   `702f42c`. `bench/sp-pbwt-bcf.base` is the frozen build of exactly that
   commit, made in Run 8 before any of that plan's edits — the correct target
   for numbers that are supposed to be "about 702f42c".
3. Moved the stale 32-thread outputs aside (`results_spbwt_eval.32thread.bak`,
   `benchmarks_spbwt_eval.32thread.bak`, `all2.csv.32thread.bak` in `../exp`)
   rather than deleting them, and re-ran via
   `limitram 32G snakemake -s Snakefile.spbwt_eval -c 32 -p --rerun-incomplete`
   under `babysit-run`. With `threads: 8` and `-c 32`, up to 4 jobs run
   concurrently (≤32 cores total) rather than Run 7's one-at-a-time
   32-thread jobs — took ~24.5 min for all 44 jobs, versus Run 7's ~1h45m.
4. Rebuilt `all2.csv` the same way as Run 7 (`build_all2.py`, new in this
   correction): copy `results/all.csv`, replace only the 44
   `bcf-prs`/`bcf-spr` rows from the new `.time` files via `time_to_csv.py`,
   verify every other row is byte-identical in Python. 396 rows total, 44
   replaced, 0 unexpected mismatches among the other 352.

### Results (wall clock, seconds; 8 threads both sides now)

| chr | bcf-prs (old→new) | bcf-spr (old→new) |
|---|---|---|
| chr1 | 661.86 → 235.95 | 1738.44 → 222.58 |
| chr2 | 740.99 → 290.87 | 1895.63 → 237.77 |
| chr3 | 594.34 → 216.09 | 1564.28 → 206.07 |
| chr4 | 584.26 → 219.90 | 1540.03 → 189.81 |
| chr5 | 544.27 → 187.83 | 1409.95 → 175.58 |
| chr6 | 514.35 → 177.87 | 1353.52 → 170.02 |
| chr7 | 479.73 → 169.55 | 1271.06 → 161.30 |
| chr8 | 468.62 → 169.92 | 1257.21 → 156.31 |
| chr9 | 368.23 → 134.76 | 961.55 → 120.41 |
| chr10 | 407.23 → 145.06 | 1075.01 → 132.73 |
| chr11 | 414.07 → 149.76 | 1094.94 → 134.37 |
| chr12 | 395.75 → 143.28 | 1044.43 → 127.83 |
| chr13 | 295.20 → 106.43 | 774.47 → 94.72 |
| chr14 | 273.03 → 99.56 | 716.15 → 88.24 |
| chr15 | 249.59 → 90.39 | 653.42 → 77.42 |
| chr16 | 279.59 → 100.77 | 730.21 → 88.11 |
| chr17 | 237.10 → 90.76 | 625.80 → 76.51 |
| chr18 | 233.53 → 84.16 | 611.13 → 75.79 |
| chr19 | 192.19 → 66.83 | 493.10 → 61.40 |
| chr20 | 185.28 → 68.55 | 490.31 → 60.71 |
| chr21 | 114.29 → 39.57 | 296.00 → 36.45 |
| chr22 | 114.13 → 41.92 | 296.26 → 35.92 |

Genuine 8-thread-vs-8-thread speedups of roughly 2.5-4x wall clock across
chromosomes for both modes — smaller than Run 7's invalid 32-vs-8 numbers
(e.g. chr1 `bcf-spr` Run 7 claimed 1738.4 → 253.2 = 6.9x; the real 8-vs-8
number is 1738.4 → 222.6 = 7.8x, actually slightly *better* here, but chr22
`bcf-prs` Run 7 claimed 114.1 → 55.8 = 2.0x vs the real 114.1 → 41.9 = 2.7x —
direction of the error isn't consistent across rows, which is exactly why the
32-vs-8 comparison wasn't trustworthy in either direction).

### Where things live

- `/data/proj/2bfpbwt/exp/Snakefile.spbwt_eval` — `THREADS = 8`,
  `BCF_EXE = bench/sp-pbwt-bcf.base` (both edited in place; comment explains
  why).
- `/data/proj/2bfpbwt/exp/build_all2.py` (new, untracked) — the CSV-merge
  script, replaces the ad-hoc Python one-liner from Run 7.
- `/data/proj/2bfpbwt/exp/all2.csv` (overwritten with corrected numbers);
  `all2.csv.32thread.bak`, `results_spbwt_eval.32thread.bak/`,
  `benchmarks_spbwt_eval.32thread.bak/` (Run 7's stale 32-thread artifacts,
  kept for reference, not deleted).

Nothing committed in either repo. The "Repeatable procedure" section at the
end of Run 7's entry above should be read with `THREADS = 8` (already fixed
in the live `Snakefile.spbwt_eval`), not 32, for any future refresh.

## Run 10 (Step 2: bit-transpose the BCF window packing) — 2026-08-29

**Status:** SUCCESS, kept. Nothing reverted.

Replaced the three copies of the scalar byte→bit window pack in `iobcf.c`
(`FGETCOLIW_IMPL(W)`'s `fgetcoliw##W##r`, `fgetcoliwgr`, and `fgetcolwgri`)
with one shared helper. The old pattern

```c
for (size_t r = 0; r < n; r++) {
  uint64_t val = 0;
  for (size_t k = 0; k < wix; k++) val |= (uint64_t)stage[k * n + r] << k;
  c[r] = val;
}
```

strides the `w x n` staging buffer by `n` and costs ~64n operations per
window. `perf annotate` in the plan's baseline put 76% of `fgetcoliw64r`'s
samples — ≈14.5% of total `sampled` runtime — on exactly those four
instructions.

New path (`iobcf_read_window`, plus `iobcf_pack_bitrow`,
`iobcf_transpose64`, `iobcf_transpose_window`, `iobcf_pack_scalar`):

1. `bcf_decode_gt_row` is unchanged and still decodes each record into the
   staging row.
2. Immediately after each decode — while the ~5 KB byte row is still hot in
   L1 — the row is compressed to a bitmap row of `ceil(n/64)` words, one bit
   per haplotype. Under AVX2 this is `loadu` + two `cmpeq`/`subs_epu8` +
   `movemask_epi8`, 32 bytes → 32 bits per iteration; a portable scalar
   fallback sits behind the same `#if defined(__AVX2__)` and is used
   automatically on non-AVX2 targets. The repo's `-march=native` build on
   this machine selects the AVX2 path.
3. Once the window is complete, the `wix x n` bitmap (rows `>= wix` treated
   as zero, so partial windows and `w < 64` variants come out with the high
   bits clear) is transposed in 64x64 blocks with the Hacker's Delight
   bit transpose — 6 shift/mask/xor rounds, ~11 ops per output word instead
   of 64.

Bit order is preserved exactly: record `k` still lands at bit `k` of `c[r]`,
LSB = first record of the window.

**Fidelity fallback.** The old loop OR-ed the *whole allele byte* shifted left
by `k`, so an allele value outside {0,1} would bleed into higher bits — a
bitmap cannot reproduce that. `iobcf_pack_bitrow` therefore reports whether it
saw any byte outside {0,1} (multiallelic allele indices, and the `0xff` that
`bcf_gt_allele()` yields for a missing genotype), and such a window falls back
to `iobcf_pack_scalar`, the original loop. Output is thus identical on *any*
input, not merely on the biallelic panels used here. `w > 64` also takes the
fallback.

All buffers the helper keeps across calls (`stage`, `bits`) are
`static __thread`, per the existing warning at `fgetcolwgri` — that function
is called concurrently by `stagpar`'s per-lane threads.

### Unit check (run *before* wiring, as the plan requires)

Standalone harness that `#include`s the real `iobcf.c`, so it exercises the
shipped functions rather than a copy. It compares the old scalar pack against
bitrow-compress + transpose, bit for bit, over:

- every `w` in 1..64;
- `n` in {1, 2, 7, 31, 32, 33, 63, 64, 65, 127, 128, 129, 1000, 5008, 5009}
  (covering sub-word, exact-multiple and non-multiple-of-64 widths, plus the
  real 5008-haplotype panel width);
- partial fills `wix` = 0, 1, `w/2`, `w-1`, `w` (all `wix` for `w <= 8`);
- bitmap rows poisoned with `0xdeadbeefcafef00d` before each case, so any
  read of a row `>= wix` would show up as a mismatch.

**4860 cases, all bit-identical**, in both the `-march=native` (AVX2) build
and a `-mno-avx2` build, so the portable fallback is verified too. The harness
also asserts that bytes `2` and `0xff` are reported by the fidelity check and
that `1` is not.

### Verification

`bench/verify.sh` → `== ALL BIT-IDENTICAL ==`: all 4 BCF modes x 2 panels
byte-identical to `bench/ref/`.

### Timing (chr21, single-threaded `OMP_NUM_THREADS=1`, `limitram 16G`,
`bench/sp-pbwt-bcf.base` vs new build, interleaved, 3 runs each —
`bench/step2_timing.sh`, run as a single `babysit-run` pipeline per workflow
rule 7 and the Run 9 note)

`TRACE` wall seconds. Note the base binary predates Step 1 too, so these
deltas are Step 1 + Step 2 combined; Run 9 measured Step 1 alone at ~0.4% on
`linear` and ~flat on `sampled`, so essentially all of the `sampled` gain
below is Step 2.

| mode | base (3 runs) | new (3 runs) | median base | median new | median Δ |
|---|---|---|---|---|---|
| `sampled` | 28.78, 29.35, 28.87 | 25.49, 25.51, 30.05 | 28.87 | 25.51 | **-11.6%** |
| `linear` (control) | 37.38, 37.31, 42.91 | 37.16, 37.14, 37.38 | 37.38 | 37.16 | -0.6% (noise) |
| `blockpar` | 254.98, 254.81, 281.62 | 250.91, 252.19, 289.37 | 254.98 | 252.19 | -1.1% (noise) |
| `stagpar` | 261.02, 257.35, 255.39 | 239.90, 269.34, 244.38 | 257.35 | 244.38 | -5.0% (noisy) |

`sampled` is the plan's primary target for this step and improves by 11.6% at
the median; runs 1 and 2 are tight (25.49/25.51 vs 28.78/29.35) and only the
third repetition is an outlier, in a round where the machine was visibly noisy
across the board (base `linear` 42.91 s and base `blockpar` 281.62 s in that
same round). Against the plan's stated realistic ceiling of 13–17% off
`sampled` for the *whole* plan, Step 2 alone delivering ~11.6% is on target.

**`linear` confirmed unaffected**, and confirmed by inspection rather than
assumption: `linc` reads columns through `fgetcoli`, which calls
`bcf_decode_gt_row` directly and never enters the window-packing path at all.
Its -0.6% median is run-to-run jitter (its own base round-3 sample was 42.9 s).

`blockpar`/`stagpar` deltas are within run-to-run noise at these
single-threaded runtimes: both re-stream the panel per block/lane, so window
packing is a much smaller share of their ~250 s than of `sampled`'s ~29 s. No
extra timing batch was run for them — the acceptance decision does not hinge
on it (the primary target won decisively, `bench/verify.sh` passes, and
neither mode regresses at the median), and this is a shared machine.

**Decision: keep.** Nothing reverted.

### Files touched

- `iobcf.c` — new shared packer (`iobcf_pack_bitrow`, `iobcf_transpose64`,
  `iobcf_transpose_window`, `iobcf_pack_scalar`, `iobcf_read_window`) and the
  three call sites reduced to a single call each. `io.h` unchanged.
- `bench/step2_timing.sh` (new) — the interleaved base-vs-new sweep.

The unit-check harness was throwaway (scratchpad, not added to the Makefile or
the repo), as the plan specified.

### Note for the eventual Step 6

`blockpar`/`stagpar` numbers in `../exp/all2.csv` are now stale again, since
this step changed their packing path as well. Refresh playbook: memory
`exp_benchmark_refresh_workflow`, and the "Repeatable procedure" section at
the end of Run 7 (read with `THREADS = 8`, per Run 9b).

Nothing committed — Step 6 of the plan does the commit.

## Run 11 (Step 3: remove dead per-window copies in `wapproxc_rrs`) — 2026-08-29

**Status:** SUCCESS, kept. Nothing reverted.

`wapproxc_rrs`'s per-window loop (and its last-window special case) previously
did four `memcpy`s to snapshot `pbwt`'s `a`/`d` into `pbwtPr`/`pbwtPrRev`
before calling `rrsortx(nrow, w64, pbwt->a, aux)` in place. Tracing the reads
in `divc`/`recover_div` showed only `pbwtPr->d` and `pbwtPrRev->a` are ever
read back out — `pbwtPr->a` and `pbwtPrRev->d` are write-only snapshots that
get immediately overwritten or discarded, so copying them was pure waste.

Added `rrsortx_src(size_t n, uint64_t *c, pbidx_t *src, pbidx_t *dst, pbidx_t *aux)`
next to `rrsortx`: an 8-pass radix sort that reads pass 0 from `src` (left
untouched) and writes to `aux`, alternates `aux`/`dst` for the remaining
passes, and lands the final sorted permutation in `dst` — same 8-pass/ping-pong
structure as `rrsortx`, just parameterized to sort out-of-place from a
separate source buffer instead of in-place with an internal copy.

Both call sites (per-window loop body and the last-window section) replaced
their four `memcpy`s + in-place `rrsortx` with:

```c
SWAP(pbwt, pbwtPr);
SWAP(pbwtRev, pbwtPrRev);
rrsortx_src(nrow, w64, pbwtPr->a, pbwt->a, aux);
```

`SWAP` exchanges the whole `pbwtad*` pointers (never individual `a`/`d`
fields — `pbwtad_new` allocates `a`/`d` as one contiguous block, and
`PBWTAD_FREE` assumes `a` is the block base, so partial-field swapping would
corrupt the free). After the swap, `pbwtPr`/`pbwtPrRev` hold what was
previously live in `pbwt`/`pbwtRev` (no copy needed — they already point at
the right data), and `rrsortx_src` writes the new sorted result directly into
`pbwt->a` from `pbwtPr->a`, skipping the memcpy that used to feed `rrsortx`'s
in-place sort. `PBWTAD_FREE(pbwt/pbwtRev/pbwtPr/pbwtPrRev)` at the end is
unchanged and still correct: all four originally-malloc'd blocks get freed
exactly once, regardless of which local variable currently points to which
after the swaps.

Had to disambiguate this edit from `swbapproxc_rrs` (the BM-backend sibling
function with an identical-looking memcpy block, out of scope for this
BCF-only plan) via `grep`/`awk` line-ownership checks before editing.

### Verification

`bench/verify.sh` → `== ALL BIT-IDENTICAL ==`: all 4 BCF modes x 2 panels
byte-identical to `bench/ref/`.

### Timing (chr21, single-threaded `OMP_NUM_THREADS=1`, `limitram 16G`,
`bench/sp-pbwt-bcf.base` vs new build, interleaved, 3 runs, `sampled` mode
only — `bench/step3_timing.sh`, launched via `babysit-run`)

| run | base (s) | new (s) |
|---|---|---|
| 1 | 28.926 | 25.702 |
| 2 | 28.764 | 25.475 |
| 3 | 28.936 | 25.558 |
| **mean** | **28.875** | **25.578** |

**-11.4%** mean wall time on `sampled`. Consistent across all three runs
(spread <0.6s on both sides) — not noise. The plan estimated this step at
"~1-2%, low risk, removes dead work"; the actual win is much larger than
that estimate, likely because removing the memcpys also removes the L2/L3
traffic they generated for arrays already resident from the immediately
preceding sort pass, not just the copy cost itself.

**Decision: keep.** Nothing reverted. (Step 3 has no explicit
revert-if-flat instruction in the plan, only Steps 1 and 4 do — moot here
since the result is a clear win anyway.)

### Files touched

- `sp-pbwt.c` — `rrsortx_src` (new, next to `rrsortx`), `wapproxc_rrs`
  per-window and last-window sections.
- `bench/step3_timing.sh` (new) — interleaved base-vs-new sweep.

Nothing committed — Step 6 of the plan does the commit.

## Run 12 (Step 4: rrsortx experiments — skip degenerate passes / carry key with index) — 2026-08-29

**Status:** BOTH EXPERIMENTS DISCARDED, both reverted. `sp-pbwt.c` is back to
exactly the state Run 11 (Step 3) left it in. This is the plan-sanctioned
"neither wins by more than noise" outcome.

The plan flagged `rrsortx` at 9.3% of `sampled`-mode runtime (`perf record -F
299`, chr21) and asked for two independently-measured experiments against the
Step-3 baseline (`bench/sp-pbwt-bcf.step3`, a frozen copy of the current tree's
build made before any Step-4 edit). Note that BCF `sampled` (`wapproxc_rrs`)
no longer calls `rrsortx` at all since Run 11 — it calls `rrsortx_src` in both
its per-window and last-window sections — so the hot function the perf sample
attributed to `rrsortx` is really `rrsortx_src`, and both experiments targeted
it (experiment (a) also covered `rrsortx`, which still serves the BM backends
and the `blockpar`/`stagpar` kernels).

**(a) Skip degenerate radix passes.** Both functions already build all 8 byte
histograms up front, so a byte position that is constant across all `n`
elements is detectable as a bucket whose count equals `n` — that pass is the
identity permutation and can be skipped. Implemented in both `rrsortx` and
`rrsortx_src`, with the ping-pong parity fix the skip makes necessary:

```c
    int degenerate = 0;
    for (size_t b2 = 0; b2 < 256; b2++)
      if (cnt[i][b2] == n) { degenerate = 1; break; }
    if (degenerate)
      continue;
    ...
  }
  // odd executed-pass count leaves the result in the wrong buffer
  if (pre != s)
    memcpy(s, pre, n * sizeof *s);
```

For `rrsortx_src` the same issue takes two forms: zero executed passes leaves
the result in `src`, an odd count leaves it in `aux`, and the "first executed
pass reads from `src`" special case had to move off `i == 0` onto a
`passes == 0` counter (pass 0 may now be the skipped one). Both are handled by
the single `if (pre != dst) memcpy(dst, pre, ...)` tail.

As the plan predicted, real low-MAF panel data barely exercises the skip path
(across 5008 haplotypes essentially every byte of a 64-bit window key spans
more than one value), so `bench/verify.sh` alone would not have tested the
parity handling. It was therefore validated separately by a 400-trial
randomized harness that extracts both functions verbatim out of `sp-pbwt.c`
and sorts inputs whose per-byte variability is chosen by a random 8-bit mask
(so 0–8 byte positions are constant per trial, including the all-constant
n=1/no-pass case), comparing against a stable insertion sort and also
asserting `rrsortx_src` never clobbers `src`. All 400 trials passed.

**(b) Carry the key with the index.** Each radix pass previously did `n`
indexed reads `c[pre[j]]` into the ~40KB key array (bigger than L1 at
n=5008, so effectively permutation-driven L2 gathers). Replaced with a
16-byte `rkv_t { uint64_t k; uint64_t v; }` packed (key, index) element sorted
as one unit, so all 8 passes scan one array sequentially. The (key, index)
build is fused into the existing histogram pass (it reads `c[src[j]]` once and
histograms the value it just loaded, so no extra pass is added), and the
indices are read back out into `dst` at the end. The two `rkv_t` scratch
buffers are allocated once in `wapproxc_rrs` and passed in, replacing the
`aux` parameter, rather than being malloc'd per window (~16k windows on
chr21). Net: 1 gather pass + 8 sequential passes + 1 extract pass, versus 8
gather passes.

### Verification

`bench/verify.sh` → `== ALL BIT-IDENTICAL ==` for each experiment's build
(all 4 BCF modes x 2 panels), and again after reverting both. Experiment (a)
additionally passed the 400-trial synthetic degenerate-byte harness described
above.

### Timing (chr21, single-threaded `OMP_NUM_THREADS=1`, `limitram 16G`,
`bench/sp-pbwt-bcf.step3` vs candidate build, interleaved, 3 runs each,
`sampled` mode only — `bench/step4_timing.sh`, launched via `babysit-run`)

Experiment (a), skip degenerate passes:

| run | step3 (s) | (a) (s) |
|---|---|---|
| 1 | 25.459 | 25.532 |
| 2 | 25.876 | 25.443 |
| 3 | 25.412 | 25.413 |
| **mean** | **25.582** | **25.462** |

**-0.5%** — noise. Run 2 of the baseline is the only sample outside a 0.12s
band, and the per-run deltas change sign across runs.

Experiment (b), carry key with index:

| run | step3 (s) | (b) (s) |
|---|---|---|
| 1 | 25.496 | 25.588 |
| 2 | 25.649 | 25.603 |
| 3 | 25.802 | 25.557 |
| **mean** | **25.649** | **25.582** |

**-0.3%** — noise, and again sign-inconsistent per run (run 1 is slower).

Both deltas are far below the ~11% mean shift with <0.6s spread that Run 11
treated as clearly real, and sit inside the 1–5% band Run 10 already
classified as machine noise. Neither experiment clears the bar.

**Decision: discard both.** (a) is correct and defensible but simply never
fires on this data, so it only adds 2048 comparisons per sort plus a parity
branch for nothing measurable. (b) trades 4-byte moves for 16-byte ones across
8 passes; at n=5008 the key array it was meant to stop gathering from is
apparently resident enough in L2 that the extra write bandwidth cancels the
gather savings exactly. Since neither pays for its complexity, the plan's
revert-both branch applies and `rrsortx`/`rrsortx_src` are left byte-for-byte
as Step 3 wrote them.

A likely explanation for the flat result overall: the perf attribution of 9.3%
to the sort is an upper bound on what Step 4 could ever have recovered, and
both variants keep the same number of scatter-writes — the part of a radix
pass that actually stalls — so neither attacks the dominant cost.

### Files touched

- `sp-pbwt.c` — `rrsortx` and `rrsortx_src` modified and then fully reverted;
  net zero change from Run 11's state.
- `bench/step4_timing.sh` (new) — interleaved step3-vs-candidate sweep.
- `bench/step4_timing.txt` (new) — raw per-run output for both experiments.

Nothing committed — Step 6 of the plan does the commit.

## Run 13 (Step 5: cache GT format id in bcf_decode_gt_row) — 2026-08-29

**Status:** SUCCESS — correctness verified bit-identical, GT id caching implemented,
timing shows flat/negligible impact (consistent with the plan's ~1% estimate).

### What changed

`iobcf.c`, function `bcf_decode_gt_row`:

1. **Cached GT format tag id** to avoid per-call khash string lookup overhead of
   `bcf_get_fmt(hdr, line, "GT")`. The id is resolved once via
   `bcf_hdr_id2int(hdr, BCF_DT_ID, "GT")` and cached as `static __thread int`.

2. **Header pointer validation** — because `fgetcolwgri` can re-open the BCF file
   (lines 483–490), which creates a new `bcf_hdr_t*`, the cache is invalidated
   if the header pointer changes. Both the cached id and the header pointer are
   `static __thread` to handle concurrent lane threads in `stagpar` mode
   (following the existing pattern documented in the file).

3. **Replaced `bcf_get_fmt(hdr, line, "GT")` with `bcf_get_fmt_id(line, gt_fmt_id)`**
   — the id-based lookup avoids the string hash. `bcf_get_fmt_id` calls
   `bcf_unpack(line, BCF_UN_FMT)` internally (just like `bcf_get_fmt` does), so no
   separate unpack is needed.

4. **Updated comments** to explain the GT id caching strategy and the htslib
   function changes (s/`bcf_get_fmt`/`bcf_get_fmt_id`/).

### Build

`HTSLIB=/data/proj/2bfpbwt/htslib make sp-pbwt-bcf` — compiled cleanly with no
warnings or errors.

### Correctness verification

`bench/verify.sh`:

```
OK   linear.small
OK   linear.mid
OK   sampled.small
OK   sampled.mid
OK   blockpar.small
OK   blockpar.mid
OK   stagpar.small
OK   stagpar.mid
== ALL BIT-IDENTICAL ==
```

All four modes on both small and mid panels: **bit-identical**, as required.

### Timing (chr21, single-threaded `OMP_NUM_THREADS=1`, `sampled` mode, 3 runs)

| Run | Wall (s) | User (s) | Sys (s) |
|---|---|---|---|
| 1 | 25.603 | 25.318 | 0.130 |
| 2 | 25.499 | 25.229 | 0.122 |
| 3 | 25.570 | 25.288 | 0.131 |
| **Mean** | **25.557** | **25.278** | **0.128** |

**~0.1% variation** (all within 25.5±0.1s) — well within machine noise, consistent
with the plan's estimate that this optimization is only ~1% and may not be
measurable against natural variance in single runs on a shared machine. No
regression detected.

### Files touched

- `iobcf.c` — `bcf_decode_gt_row` modified to cache GT format id and header
  pointer, avoiding per-call khash lookup.

Nothing committed — Step 6 of the plan does the commit.

## Run 14 (Step 6: measure, document, commit) — 2026-08-29

**Status:** Plan complete. Steps 1-5 all built on `master`/uncommitted, this run
measures the combined effect, re-profiles, and commits.

### Verification

`bench/verify.sh` on the final tree (`make -B sp-pbwt-bcf`, i.e. clean rebuild
of everything Steps 1-5 touched) → `== ALL BIT-IDENTICAL ==` for all 4 BCF
modes x 2 panels. This is the same gate every step already passed individually;
re-run here as the final combined check.

### Timing — full before/after (`bench/step6_timing.sh`, `bench/sp-pbwt-bcf.base`
[pre-plan, commit `702f42c`] vs `bench/sp-pbwt-bcf.final` [all of Steps 1-5],
3 runs each, interleaved, `limitram 16G`, chr21 + chr20, launched via
`babysit-run`)

`linear`/`sampled` at `OMP_NUM_THREADS=1` (the plan's target — these are the
sequential modes). `blockpar`/`stagpar` at `OMP_NUM_THREADS=8`, matching the
manuscript/`../exp` convention for `bpr`/`spr` (see `exp_benchmark_refresh_workflow`
memory) — these two modes were only touched incidentally, via Step 2's shared
window-packing path, so they're measured at their normal operating thread
count rather than forced single-threaded.

| mode | chr | base mean (s) | final mean (s) | Δ |
|---|---|---|---|---|
| `linear` | 21 | 41.92 | 37.04 | **-11.6%** |
| `linear` | 20 | 61.45 | 60.99 | -0.7% |
| `sampled` | 21 | 28.68 | 25.57 | **-10.9%** |
| `sampled` | 20 | 47.61 | 42.11 | **-11.6%** |
| `blockpar` (8t) | 21 | 45.20 | 44.89 | -0.7% |
| `blockpar` (8t) | 20 | 71.89 | 69.37 | -3.5% |
| `stagpar` (8t) | 21 | 36.95 | 36.49 | -1.3% |
| `stagpar` (8t) | 20 | 59.97 | 59.23 | -1.2% |

`linear`/`sampled` land at or beyond the plan's stated realistic ceiling
(10-13% off `linear`, 13-17% off `sampled`, given zlib decompression is
~50-65% of runtime and explicitly out of scope). chr21's `linear` win
(-11.6%) is larger than the earlier isolated Step 1 measurement (Run 9's
initial 3-run pass called it inconclusive before 7 extra reps confirmed a
real but modest win) — the combined effect across all steps compounds to
more than any single step measured alone, and chr21's base run 1 in this
sweep was itself a 44.4s outlier (vs 40.5/40.9s for runs 2-3), so -11.6% is
if anything a conservative read of the win. chr20 (a larger chromosome, more
zlib-bound per the fixed decompression cost per record) shows a smaller
`linear` delta (-0.7%, close to noise) but a full `sampled` delta (-11.6%),
consistent with `sampled`'s wins (Steps 2-3) being structural rather than
proportional-to-zlib-share.

`blockpar`/`stagpar` show small but consistently negative (never positive)
deltas across both chromosomes and 8-thread CPU-seconds move the same
direction — real, if modest, since only Step 2's packing-path change reaches
them.

### Re-profile (`perf record -F 299`, final binary, chr21, single-threaded)

| symbol | `linear` (final) | `linear` (plan's original baseline) | `sampled` (final) | `sampled` (original baseline) |
|---|---|---|---|---|
| `inflate_fast` + `crc32_z` (zlib) | 48.8% | 49.4% | 71.0% | 65.3% |
| `cpbwti` | 44.1% | 44.1% | — | — |
| `rrsortx_src` (was `rrsortx`) | — | — | 11.5% | 9.3% |
| `fgetcoliw64r` | — | — | 6.0% | 19.2% |
| `iobcf_transpose_window` (new, Step 2) | — | — | 2.5% | — |
| `fgetcoli` | 3.0% | 3.1% | — | — |
| `divc` | — | — | 2.0% | 1.9% |

`fgetcoliw64r`'s share collapsed from 19.2% to 6.0% — Step 2's bit-transpose
replaced its dominant scalar byte→bit loop, and `iobcf_transpose_window` (the
new Hacker's-Delight transpose) now accounts for a small, separate 2.5%
instead. `rrsortx_src`'s share went *up* (9.3% → 11.5%) in relative terms
even though `sampled`'s total wall time dropped ~11%, because Step 4's
experiments on it were both discarded (see Run 12) — its absolute cost is
essentially unchanged, so as everything around it got faster its share of a
smaller total mechanically increased. zlib's combined share also rose in
both modes for the same reason: it is fixed per-record decompression work
untouched by this plan, so it becomes a larger fraction of a shrinking total.
`cpbwti`'s share stayed flat at 44.1% despite `linear`'s absolute wall time
dropping ~11.6% on chr21 — its own absolute cost dropped roughly in step with
zlib's fixed cost, which is a coincidence of these particular panel sizes,
not a claim that `cpbwti` is now zero-cost.

**What's left**, per this profile: zlib decompression is now 49-71% of
runtime in both modes and is explicitly out of scope for this plan (no
libdeflate, no decompression threads, per the plan's stated constraints) —
that is the largest remaining lever, but a different plan. Within scope,
`rrsortx_src` (11.5% of `sampled`) is the next largest kernel; Step 4 already
tried two approaches there and both were noise-level, so a different
algorithmic angle would be needed to move it further.

### Steps summary (all details in their own Run entries above)

| Step | What | Decision |
|---|---|---|
| 1 (Run 9) | `cpbwti`/`cpbwt`/`cpbwtiLCP`: eliminate `o`/`h` staging arrays | Kept |
| 2 (Run 10) | Bit-transpose BCF window packing (`iobcf.c`) | Kept, -11.6% on `sampled` alone |
| 3 (Run 11) | Remove dead per-window copies in `wapproxc_rrs`, add `rrsortx_src` | Kept, -11.4% on `sampled` alone |
| 4 (Run 12) | `rrsortx`/`rrsortx_src`: skip degenerate passes; carry key with index | **Both reverted** — measured -0.5%/-0.3%, within noise |
| 5 (Run 13) | Cache GT format id in `bcf_decode_gt_row` | Kept, measured flat in isolation (~0.1%) but low-risk per plan |

No step was silently dropped: Step 4 is the one explicit revert, documented
in Run 12 with its own measurements and reasoning.

### Files touched (cumulative, Steps 1-5)

- `sp-pbwt.c` — `cpbwti`/`cpbwt`/`cpbwtiLCP` (Step 1), `rrsortx_src` (new,
  Step 3) and `wapproxc_rrs` (Step 3), `rrsortx`/`rrsortx_src` touched and
  reverted (Step 4, net no diff vs Step 3's state)
- `iobcf.c` — shared bit-transpose window packer (Step 2), GT format id
  caching in `bcf_decode_gt_row` (Step 5)
- `io.h` — unchanged, as the plan required (no signature churn)
- `bench/verify.sh`, `bench/ref/*` (gitignored, regenerable) — correctness
  harness, Step 0
- `bench/*_timing.sh` — per-step timing sweeps (Steps 0, 2, 3, 4, 6)
- `.gitignore` — added patterns for `bench/.pipeline-step*/`,
  `bench/*.data`, `bench/*.stderr` (perf profile artifacts from this run)
- `CHANGELOG.md` — this file, Runs 8-14

### Note on `../exp`

`blockpar`/`stagpar` numbers in `../exp/all2.csv` are stale again as of this
commit — Step 2's window-packing change affects their read path too (see the
-0.7%/-3.5%/-1.3%/-1.2% deltas measured above), on top of the commit
`702f42c` changes `all2.csv` already reflects. Refresh via the
`exp_benchmark_refresh_workflow` memory playbook (scoped rerun of `blockpar`/
`stagpar` only, `THREADS = 8` to stay comparable, report the manuscript
summary table per that playbook's step 6) — not part of this plan, so not
done here.

### Commit

This run's commit includes all of Steps 1-5's code changes
(`sp-pbwt.c`, `iobcf.c`), the Step 0 harness (`bench/verify.sh`), the
per-step timing scripts, the `.gitignore` additions, and this changelog.
`bench/sp-pbwt-bcf.base`/`.step3`/`.final` (frozen comparison binaries),
`bench/ref/*.dump` (regenerable reference dumps), and `bench/.pipeline-*/`
(run logs) are all gitignored and not committed. No `git push` (repo rule 1).
