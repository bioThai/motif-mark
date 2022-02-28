#!/bin/bash

#need this before activating conda env if running shell script locally
eval "$(conda shell.bash hook)"

conda activate motif_mark

input_fasta_dir="./test_input_files/"
input_motif_dir="./test_input_files/"
input_fasta_file="Figure_1.fasta"
input_motif_file="Fig_1_motifs.txt"

python3 motif-mark-oop.py -f $input_fasta_dir${input_fasta_file} -m $input_motif_dir${input_motif_file}


exit