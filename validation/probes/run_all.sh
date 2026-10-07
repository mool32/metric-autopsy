#!/usr/bin/env bash
# Run every exploratory probe against the checkout's engine (~3 min).
# p12 (memory) is opt-in because it peaks at ~6 GB RAM:  WITH_MEMORY=1 ./run_all.sh
set -e
cd "$(dirname "$0")"
echo "# metric-autopsy probes — git $(git rev-parse --short HEAD 2>/dev/null || echo '?')," \
     "python $(python -c 'import platform; print(platform.python_version())')," \
     "numpy $(python -c 'import numpy; print(numpy.__version__)')"
for p in p0*.py p1[01]*.py; do
    echo "=================== $p"
    python "$p"
done
if [ "${WITH_MEMORY:-0}" = "1" ]; then
    echo "=================== p12_memory_scale.py"
    python p12_memory_scale.py 5000
fi
