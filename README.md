# Validation v1: results of workflow run 38070771737

The frozen code is the tag v0.3.0-prereg; the run tag panel-v1-run names the drand round whose randomness is the key (key.json).
reports.tar.gz holds one JSON report per claim card; SHA256SUMS covers every file of the uncompressed results.
compact.tar (or its parts, concatenated) holds the planned backgrounds the shards read (compact/SHA256 is their digest); pilot.json and backgrounds.json are the run tag's.
Anyone can re-run datasets and compare the reports: tar -xf compact.tar and cat reports.tar.gz* | tar -xzf - (concatenate the parts first where they are split), then, at the run tag, python validation/prereg/blind.py verify --results . --compact compact --pilot pilot.json --datasets 20 --out rerun (it checks SHA256SUMS first)
