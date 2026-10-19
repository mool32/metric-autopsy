import sys
# usage: add.py MODE(w|a) < lines with fields separated by " | "
OUT = "/tmp/claude-0/-home-user-metric-autopsy/db456a51-1a6c-55a4-a0ce-5b88600b19d7/scratchpad/lb_screen/p2_166ad225/batch_15.tsv"
mode = sys.argv[1]
rows = []
for line in sys.stdin.read().splitlines():
    if not line.strip():
        continue
    parts = [p.strip() for p in line.split(" | ")]
    assert len(parts) == 5, line
    rid, dec, code, kinds, note = parts
    assert dec in ("KEEP", "EXCLUDE", "UNCLEAR"), line
    if dec == "EXCLUDE":
        assert code in ["X%d" % i for i in range(1, 11)], line
    else:
        assert code == "-", line
    assert kinds == "-" or all(k in ("pair", "set", "global") for k in kinds.split(",")), line
    assert "\t" not in note and len(note.split()) <= 30, (rid, len(note.split()))
    rows.append((rid, dec, code, kinds, note))
with open(OUT, mode) as f:
    if mode == "w":
        f.write("rid\tdecision\tcode\tkinds\tnote\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
print("wrote", len(rows), "last", rows[-1][0] if rows else None)
