# SPDX-License-Identifier: MIT
"""Compare regenerated result JSON files with the shipped ones, across platforms.

Usage: python -I -B tools/ci_compare.py NEW_DIR REFERENCE_DIR NAME.json [NAME.json ...]

The `_environment` block is ignored. Every other value must be equal, except that
floating-point numbers may differ by a relative tolerance of 1e-6 (absolute 1e-12),
which absorbs last-digit differences between platforms' maths libraries. Exact
byte-identity of the JSON (apart from `_environment`) is also reported, because
the macOS release replay is expected to be identical. Exit 1 on any difference
beyond tolerance, on a missing file, or on a non-empty `_failures` list.
"""
import json
import math
import sys

REL, ABS = 1e-6, 1e-12


def walk(a, b, path, bad):
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            bad.append(f"{path}: keys differ {sorted(set(a) ^ set(b))}")
            return
        for k in a:
            walk(a[k], b[k], f"{path}.{k}", bad)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            bad.append(f"{path}: lengths {len(a)} != {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", bad)
    elif isinstance(a, bool) or isinstance(b, bool) or not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        if a != b:
            bad.append(f"{path}: {a!r} != {b!r}")
    else:
        if not (math.isclose(a, b, rel_tol=REL, abs_tol=ABS) or (math.isnan(a) and math.isnan(b))):
            bad.append(f"{path}: {a!r} != {b!r}")


def main():
    new, ref, names = sys.argv[1], sys.argv[2], sys.argv[3:]
    ok = True
    for name in names:
        try:
            a = json.load(open(f"{new}/{name}"))
            b = json.load(open(f"{ref}/{name}"))
        except OSError as exc:
            print(f"MISSING {name}: {exc}")
            ok = False
            continue
        a.pop("_environment", None)
        b.pop("_environment", None)
        bad = []
        if a.get("_failures"):
            bad.append(f"_failures not empty: {a['_failures']}")
        walk(a, b, name, bad)
        exact = json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
        if bad:
            ok = False
            print(f"DIFFERENT {name}: {len(bad)} value(s) beyond tolerance")
            for line in bad[:20]:
                print("   ", line)
        else:
            print(("IDENTICAL " if exact else "WITHIN-TOLERANCE ") + name)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
