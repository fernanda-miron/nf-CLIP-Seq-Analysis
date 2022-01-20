#!/bin/bash

rm -rf test/results/  #borra resultado de prueba anterior
mkdir -p test/results/

# Correr el modulo
bash runmk.sh \
&& mv ./test/data/*.sorted.bam \
  test/results \
&& echo "=== Prueba de modulo exitosa ==="