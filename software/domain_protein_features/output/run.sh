#!/bin/bash

mkdir output
echo "Time: $(date)." >> output/info.log
echo "python scripts/run.py -hp 1 -n 1 -o output"  >> output/info.log
python scripts/run.py -hp 1 -n 1 -o output  >> output/info.log