"""Compare replay JSON outputs with the reference outputs, ignoring the _environment block."""
import json, sys
ok = True
for name in sys.argv[3:]:
    a = json.load(open(f"{sys.argv[1]}/{name}")); b = json.load(open(f"{sys.argv[2]}/{name}"))
    a.pop("_environment", None); b.pop("_environment", None)
    same = json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    print(("IDENTICAL " if same else "DIFFERENT ") + name)
    ok &= same
sys.exit(0 if ok else 1)
