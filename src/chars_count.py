#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from Bio import SeqIO
import sys

def chrom_intervals_dict(filename):
    chrom_intervals = {}
    with open(filename) as f:
        for line in f:
            name, start, stop = line.split()
            chrom_interval = chrom_intervals.get(name, [])
            chrom_interval.append(tuple([int(start), int(stop)]))
            chrom_intervals[name] = chrom_interval        
    return chrom_intervals
    
def calc_counts(str, chars_names):
    str = str.upper()
    dict ={}
    dict = dict.fromkeys(chars_names, 0)
    for i in xrange(0, len(str)-2):
        nuc_name = str[i] # A, C, T, G
        two_mer_name = nuc_name + str[i+1] # AA, AC, .., GT, GG
        three_mer_name = two_mer_name + str[i+2] # AAA, AAC, ..., GGT, GGG

        dict[nuc_name] = dict.get(nuc_name, 0) + 1;
        dict[two_mer_name] = dict.get(two_mer_name, 0) + 1;
        dict[three_mer_name] = dict.get(three_mer_name, 0) + 1;

    prevlast_nuc_name = str[len(str)-2]
    last_two_mer_name = prevlast_nuc_name + str[len(str)-1]
    last_nuc_name = str[len(str)-1]
    dict[prevlast_nuc_name] = dict.get(prevlast_nuc_name, 0) + 1;
    dict[last_two_mer_name] = dict.get(last_two_mer_name, 0) + 1;
    dict[last_nuc_name] = dict.get(last_nuc_name, 0) + 1;
    # Добавляем АТ и GC:
    dict['A/T'] = dict.get('A', 0) + dict.get('T', 0)
    dict['G/C'] = dict.get('G', 0) + dict.get('C', 0)
    dict['A/T/G'] = dict.get('A', 0) + dict.get('T', 0) + dict.get('G', 0)
    dict['T/G/C'] = dict.get('T', 0) + dict.get('G', 0) + dict.get('C', 0)
    return dict


def get_charachteristics_names():
    names = ['A', 'C', 'T', 'G']
    chars = []
    for n in names:
        chars.append(n)
        for n2 in names:
            chars.append(n + n2)
            for n3 in names:
                chars.append(n + n2 + n3)
    sortedKmers = sorted(chars, key = lambda k: (len(k), k))
    sortedKmers.insert(4, 'G/C')
    sortedKmers.insert(5, 'A/T')
    sortedKmers.insert(6, 'A/T/G')
    sortedKmers.insert(7, 'T/G/C')
    return sortedKmers

fasta_file = sys.argv[1]
bed_file = sys.argv[2]
out_file = sys.argv[3]
out_file2 = sys.argv[4]
chroms_intervals = chrom_intervals_dict(bed_file)
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
chrom_name_def = "Chromosome: "
with open(out_file, 'w') as out, open(out_file2, 'w') as out2:
    # Записываем названия колонок для итогового файла
    charachteristics_names = get_charachteristics_names()
    columns = str(charachteristics_names).replace('[', '')
    columns = columns.replace(']', '')
    columns = columns.replace('\'', '')
    columns = columns.replace(', ', '\t')
    out.write('chrom\tstart of bin\tend of bin\t' + columns + '\n')
    out2.write('nucleotides\tfrequency\tbins\t\n')
    
    for fasta in fasta_sequences:
        # Достаем имя хромосомы и последовательность
        chrom = fasta.description
        #index = chrom.index(chrom_name_def) + len(chrom_name_def)
        name = fasta.name #chrom[index : index + chrom[index:].index(" ")]
        sequence = str(fasta.seq)
        
        for chrom_interval in chroms_intervals[name]:
            start = chrom_interval[0]
            stop = chrom_interval[1]
            charachteristics = calc_counts(sequence[start : stop], charachteristics_names) # Считаем количество
            length = (stop-start) - sequence[start : stop].count('n') - sequence[start : stop].count('N') + 1
            charachteristics = {k: round((v * 1.0) / (length-len(k)), 2) for k, v in charachteristics.iteritems()} # Считаем в %
            
            for (key, value) in charachteristics.iteritems():
                row2 = key + '\t' + str(value) + '\t' + str(start) + '_' + str(stop) + '\t\n'
                out2.write(row2)

            row = name + '\t' + str(start) + '\t' + str(stop)
            for char_name in charachteristics_names:
                row = row + '\t' + str(charachteristics.get(char_name, 0))
            out.write(row + '\n')
            