#!/usr/bin/env bash
set -euo pipefail

if (( $# != 1 )); then
    echo "usage: $0 INPUT.csv" >&2
    exit 2
fi

sed -i \
    -e 's/,lrv,/,PSD+AGG,/g' \
    -e 's/,pcon,/,PSD-AGG,/g' \
    -e 's/,con,/,HC,/g' \
    "$1"
