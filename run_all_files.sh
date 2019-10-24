#!/bin/bash

echo "Reading benchmarks from $1"
echo "Looking for benchmarks in $2"
mkdir -p "$3"
while IFS='' read -r filename || [[ -n "$filename" ]]; do
    python3 main.py "$filename" "$2" "$3" >> "$4"
    echo $filename
done < "$1"
echo "Connectivity compliant circuits generated to $3"
echo "Summary file at $4"

