import sys

CHECKA = True
CHECKD = False


def file_ix(f, fname=""):
    d = dict()
    tell = f.tell()
    line = f.readline().strip()
    _c = 0
    while len(line) > 0:
        print(f"Indexing {_c}", end="\r")
        d[int(line.split(":")[0])] = tell
        tell = f.tell()
        line = f.readline().strip()
        _c += 1
    print(f"Indexing {fname} done.")
    return d


def main(_f1: str, _f2: str):
    print(_f1, _f2)
    checked = 0
    checked_wrong = 0
    f1 = open(_f1)
    idx1 = file_ix(f1, _f1)
    f2 = open(_f2)
    idx2 = file_ix(f2, _f2)

    if len(idx1) < len(idx2):
        sm = idx1
        smf = f1
        lg = idx2
        lgf = f2
    else:
        sm = idx2
        smf = f2
        lg = idx1
        lgf = f1

    for ix in sm:
        print(f"Checking {checked + 1}/{len(sm)}", end="\r")
        smf.seek(sm[ix])
        l1 = smf.readline().strip()
        lgf.seek(lg[ix])
        l2 = lgf.readline().strip()
        ix1, _data1 = l1.split(":")
        ix2, _data2 = l2.split(":")
        assert ix1 == ix2 == str(ix)

        a1, d1 = _data1.split("|")
        a2, d2 = _data2.split("|")

        _wrong = False
        msg = ""
        if CHECKA and a1 != a2:
            _wrong = True
            msg += f"{ix1}\n"
            msg += f"\tA: {a1} != {a2}\n"
        if CHECKD and d1 != d2:
            _wrong = True
            if not msg:
                msg += f"{ix1}\n"
            msg += f"\tD: {d1} != {d2}\n"
        print(msg, end="", file=sys.stderr)
        checked_wrong += 1 if _wrong else 0

        checked += 1
    print("\nDone.")
    print(f"checked: {checked}")
    print(f"  wrong: {checked_wrong}")
    print(f"  lg-sm: {set(lg) - set(sm)}")


if __name__ == "__main__":
    _f1 = sys.argv[1]
    _f2 = sys.argv[2]
    main(_f1, _f2)
