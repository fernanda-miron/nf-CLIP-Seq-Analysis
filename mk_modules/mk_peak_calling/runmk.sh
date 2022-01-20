#!/usr/bin/env bash

find -L . \
  -type f \
  -name "*.bam" \
| sed 's#.bam#.bed#' \
| xargs mk
