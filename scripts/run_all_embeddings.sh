#!/bin/bash
set -e
# override with PY=/path/to/python if python is not on your PATH
PY="${PY:-python}"
cd "$(dirname "$0")"
mkdir -p ../logs
for model in clinicalbert sapbert mpnet; do
  for clean in base stripped; do
    echo "=== $model / $clean ==="
    "$PY" generate_embeddings.py --model "$model" --clean "$clean" > "../logs/embed_${model}_${clean}.log" 2>&1
    echo "=== done $model / $clean ==="
  done
done
echo "ALL EMBEDDINGS DONE"
