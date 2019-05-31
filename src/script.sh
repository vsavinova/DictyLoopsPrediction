#!/bin/bash
python genome_size_maker.py $1 ../data/background/sizes.txt

bedtools makewindows -g ../data/background/sizes.txt -w 2000  > ../data/background/windows.bed

python chars_count.py $1 ../data/background/windows.bed $2 $3