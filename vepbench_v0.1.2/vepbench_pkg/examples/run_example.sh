#!/bin/sh
# Runs the full evaluation on the bundled simulated data (no real biological data).
python -m vepbench evaluate \
  --input examples/simulated_benchmark.csv \
  --group gene --label label \
  --scores within_signal,gene_prior,mixed,random \
  --external examples/simulated_external_ranking.csv \
  --n-boot 200 --seed 0 --out example_output
