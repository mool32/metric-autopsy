#!/usr/bin/env bash
# The inputs of the diagnostics (validation/exploratory/real_pilot_diagnostics_plan.md): the second run's
# files from results/real-pilot-2 at bad785a (case1.npz from its parts), each checked against that run's
# SHA256SUMS, and STRING v12.0's two mouse files, downloaded and tested. Nothing is computed on the data.
#
#   bash validation/exploratory/real_pilot_3_inputs.sh OUT_DIR
set -euo pipefail
out="$1"
run2=bad785a3c56e9de1ee24737953c1b4b63112aa7c
ref=refs/remotes/origin/results/real-pilot-2
git fetch -q origin "refs/heads/results/real-pilot-2:$ref"
test "$(git rev-parse "$ref")" = "$run2" || { echo "results/real-pilot-2 is not at $run2"; exit 1; }
mkdir -p "$out/data" "$out/string"
blob() { git cat-file blob "$run2:$1"; }
blob SHA256SUMS > "$out/run2-SHA256SUMS"
blob select/select.json > "$out/data/select.json"
blob run/real_pilot.json > "$out/data/real_pilot.json"
blob data/case2.npz > "$out/data/case2.npz"
: > "$out/data/case1.npz"
for p in $(git ls-tree --name-only "$run2" data/ | grep '^data/case1\.npz\.part-' | LC_ALL=C sort); do
  blob "$p" >> "$out/data/case1.npz"
done
for f in case1.npz case2.npz select.json real_pilot.json; do
  want=$(grep -E "  \./(data|select|run)/${f//./\\.}\$" "$out/run2-SHA256SUMS" | cut -d' ' -f1)
  got=$(sha256sum "$out/data/$f" | cut -d' ' -f1)
  if [ -z "$want" ] || [ "$got" != "$want" ]; then echo "$f: sha256 $got, the second run's SHA256SUMS has ${want:-nothing}"; exit 1; fi
  echo "$f: sha256 $got, as the second run's SHA256SUMS"
done
base=https://stringdb-downloads.org/download
for f in protein.info.v12.0/10090.protein.info.v12.0.txt.gz protein.links.v12.0/10090.protein.links.v12.0.txt.gz; do
  curl -fsSL --retry 5 --retry-delay 15 -o "$out/string/$(basename "$f")" "$base/$f"
  gzip -t "$out/string/$(basename "$f")"
done
(cd "$out/string" && sha256sum ./*.gz)
