#!/usr/bin/env bash

find -L . \
  -type f \
  -name "*.fastq.gz" \
| sed 's#.fastq.gz#.sam#' \
| xargs mk
