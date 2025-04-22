import sys

CHECKA = True
CHECKD = False


def main(_f1: str, _f2: str):
    print(_f1, _f2)
    checked = 0
    checked_real = 0
    checked_wrong = 0
    with open(_f1) as f1:
        with open(_f2) as f2:
            l1 = f1.readline().strip()
            l2 = f2.readline().strip()

            # print("1",l1)
            # print("2",l2)
            while l1 and l2:
                # print(l1, l2)
                _wrong = False
                ix1, _data1 = l1.split(":")
                ix2, _data2 = l2.split(":")
                assert ix1 == ix2
                if _data1 != "NULL" and _data2 != "NULL":
                    a1, d1 = _data1.split("|")
                    a2, d2 = _data2.split("|")
                    msg = ""
                    if CHECKA and a1 != a2:
                        _wrong=True
                        msg += f"{ix1}\n"
                        msg += f"\tA: {a1} != {a2}\n"
                    if CHECKD and d1 != d2:
                        _wrong=True
                        if not msg:
                            msg += f"{ix1}\n"
                        msg += f"\tD: {d1} != {d2}\n"
                    print(msg, end="", file=sys.stderr)
                    checked_real += 1
                    checked_wrong += 1 if _wrong else 0

                checked += 1
                l1 = f1.readline().strip()
                l2 = f2.readline().strip()
    print(f"checked: {checked} ({checked_real})")
    print(f"  wrong: {checked_wrong}")


if __name__ == "__main__":
    _f1 = sys.argv[1]
    _f2 = sys.argv[2]
    main(_f1, _f2)
