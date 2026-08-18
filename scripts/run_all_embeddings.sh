#!/bin/bash
set -e
PY="/c/Users/IRFANM/Downloads/pyembed/python.exe"
cd "$(dirname "$0")"
for model in clinicalbert sapbert mpnet; do
  for clean in base stripped; do
    echo "=== $model / $clean ==="
    "$PY" generate_embeddings.py --model "$model" --clean "$clean" > "../logs/embed_${model}_${clean}.log" 2>&1
    echo "=== done $model / $clean ==="
  done
done
echo "ALL EMBEDDINGS DONE"
