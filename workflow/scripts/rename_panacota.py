#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Rename sequences in a panacota file"""

from Bio import SeqIO
import glob

dic = {}
with open(snakemake.input.listfile) as f:
    for line in f.readlines()[1:]:
        name = line.split()[0]
        replacement = line.split()[1].replace(".fasta","")
        dic[name] = replacement


for file in glob.glob(f"{snakemake.input.dir_al}/*.gen"):
    fam = list(SeqIO.parse(file, "fasta"))
    for seq in fam:
        seq.id = dic[".".join(seq.id.split(".")[:3])]
        seq.description = ""
    SeqIO.write(fam, file, "fasta")

f = open(snakemake.output.done, "w")
f.write("done")
f.close()