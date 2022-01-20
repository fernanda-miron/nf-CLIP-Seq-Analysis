#!/usr/bin/env bash

find -L . \
  -type f \
  -name "*.sam" \
| sed 's#.sam#.bam#' \
| xargs mk
