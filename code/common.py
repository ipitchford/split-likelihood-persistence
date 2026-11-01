"""Shared helpers for the split-likelihood persistence certificates.

Checks use check(), never assert, so they stay active under python -O.
"""
import json
import sys
import hashlib
import platform

FAILURES = []


def check(cond, label, detail=""):
    """Record a pass/fail; failures are collected and cause a non-zero exit."""
    ok = bool(cond)
    print(("PASS " if ok else "FAIL ") + label + ((" | " + str(detail)) if detail else ""))
    if not ok:
        FAILURES.append(label)
    return ok


def expect_failure(cond_that_should_be_false, label, detail=""):
    """Negative control: the wrong claim must be detected as false."""
    return check(not bool(cond_that_should_be_false), "NEGATIVE CONTROL detects: " + label, detail)


def finish(results, path):
    """Write results JSON and exit non-zero if any check failed."""
    results["_environment"] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "optimize_flag": sys.flags.optimize,
    }
    results["_failures"] = list(FAILURES)
    with open(path, "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=True, default=float)
    with open(path, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    print("wrote", path, "sha256", digest)
    if FAILURES:
        print("FAILED CHECKS:", len(FAILURES))
        sys.exit(1)
    print("ALL CHECKS PASSED")
