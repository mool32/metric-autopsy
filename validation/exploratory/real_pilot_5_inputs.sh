#!/usr/bin/env bash
# The inputs of the third diagnostics (validation/exploratory/real_pilot_diagnostics_3_plan.md):
# - the second run's case files (results/real-pilot-2 at bad785a: case1.npz from its parts, case2.npz,
#   select.json) and the first run's case 3 (results/real-pilot at 6133dfd: case3.npz, its select.json as
#   select1.json), each checked against its run's SHA256SUMS;
# - Marsh et al. 2022's GEO records from results/marsh-metadata at 757c2ab (GSE152183's SOFT and file
#   list);
# - from GEO: GSE152183_RAW.tar (its size against the file list), and the SOFT records of GSE132042 and
#   GSE109774 (Tabula Muris Senis, Tabula Muris) for check (b4).
# Nothing is computed on the data.
#
#   bash validation/exploratory/real_pilot_5_inputs.sh OUT_DIR [--no-geo-data]
set -euo pipefail
out="$1"
geo_data=1
[ "${2:-}" = "--no-geo-data" ] && geo_data=0
run1=6133dfdf4c05f083f012aaec4136e376715af05a
run2=bad785a3c56e9de1ee24737953c1b4b63112aa7c
meta=757c2ab48af792f63b1992c8eb98b05a17d4272f
for spec in "results/real-pilot:$run1" "results/real-pilot-2:$run2" "results/marsh-metadata:$meta"; do
  br="${spec%%:*}"; want="${spec##*:}"
  git fetch -q origin "refs/heads/$br:refs/remotes/origin/$br"
  test "$(git rev-parse "refs/remotes/origin/$br")" = "$want" || { echo "$br is not at $want"; exit 1; }
done
mkdir -p "$out/data" "$out/geo"
blob() { git cat-file blob "$1:$2"; }
blob "$run2" SHA256SUMS > "$out/run2-SHA256SUMS"
blob "$run1" SHA256SUMS > "$out/run1-SHA256SUMS"
blob "$run2" select/select.json > "$out/data/select.json"
blob "$run2" data/case2.npz > "$out/data/case2.npz"
: > "$out/data/case1.npz"
for p in $(git ls-tree --name-only "$run2" data/ | grep '^data/case1\.npz\.part-' | LC_ALL=C sort); do
  blob "$run2" "$p" >> "$out/data/case1.npz"
done
blob "$run1" select/select.json > "$out/data/select1.json"
blob "$run1" data/case3.npz > "$out/data/case3.npz"
check() { # file in data/, the run's sums, the path the sums give it
  want=$(grep -E "  \./$3\$" "$2" | cut -d' ' -f1)
  got=$(sha256sum "$out/data/$1" | cut -d' ' -f1)
  if [ -z "$want" ] || [ "$got" != "$want" ]; then echo "$1: sha256 $got, its run's SHA256SUMS has ${want:-nothing}"; exit 1; fi
  echo "$1: sha256 $got, as its run's SHA256SUMS"
}
check case1.npz "$out/run2-SHA256SUMS" "data/case1\.npz"
check case2.npz "$out/run2-SHA256SUMS" "data/case2\.npz"
check select.json "$out/run2-SHA256SUMS" "select/select\.json"
check case3.npz "$out/run1-SHA256SUMS" "data/case3\.npz"
check select1.json "$out/run1-SHA256SUMS" "select/select\.json"
blob "$meta" geo/GSE152183_family.soft.gz > "$out/geo/GSE152183_family.soft.gz"
blob "$meta" geo/GSE152183_suppl_filelist.txt > "$out/geo/GSE152183_suppl_filelist.txt"
(cd "$out/geo" && sha256sum GSE152183_family.soft.gz GSE152183_suppl_filelist.txt)
get() { # url, file
  for i in 1 2 3 4; do
    if curl -fsSL --retry 3 --retry-delay 10 -m 1800 -o "$2" "$1"; then echo "ok $(stat -c %s "$2") bytes $1"; return 0; fi
    sleep $((i * 15))
  done
  echo "FAIL $1"; return 1
}
base=https://ftp.ncbi.nlm.nih.gov/geo/series
get "$base/GSE132nnn/GSE132042/soft/GSE132042_family.soft.gz" "$out/geo/GSE132042_family.soft.gz"
get "$base/GSE109nnn/GSE109774/soft/GSE109774_family.soft.gz" "$out/geo/GSE109774_family.soft.gz"
gzip -t "$out/geo/GSE132042_family.soft.gz" "$out/geo/GSE109774_family.soft.gz"
if [ "$geo_data" = 1 ]; then
  get "$base/GSE152nnn/GSE152183/suppl/GSE152183_RAW.tar" "$out/geo/GSE152183_RAW.tar"
  want=$(awk -F'\t' '$2 == "GSE152183_RAW.tar" {print $4}' "$out/geo/GSE152183_suppl_filelist.txt")
  got=$(stat -c %s "$out/geo/GSE152183_RAW.tar")
  test "$got" = "$want" || { echo "GSE152183_RAW.tar: $got bytes, GEO's file list has $want"; exit 1; }
  echo "GSE152183_RAW.tar: $got bytes, as GEO's file list"
fi
(cd "$out/geo" && sha256sum ./*)
