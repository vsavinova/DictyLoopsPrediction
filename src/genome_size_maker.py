#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from Bio import SeqIO
import sys

fasta_file = sys.argv[1]
sizes_file = sys.argv[2]
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(sizes_file, "w") as chrom_sizes_file:
    chrom_sizes = {}
    for fasta in fasta_sequences:
        # Достаем имя хромосомы и последовательность
        chrom_name = fasta.name
        chrom_seq = fasta.seq
        # Достаем размер
        chrom_size = len(chrom_seq)
        chrom_sizes[chrom_name] = chrom_size
        # Записываем в файл
        chrom_sizes_file.write(chrom_name + "\t" + str(chrom_size) + "\n")