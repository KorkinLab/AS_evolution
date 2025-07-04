#!/bin/bash

gzip -d -k -f -q *.gz
gff3_file_name=$(basename *.gff3)
fasta_file_name=$(basename *.fa)
python3 ../parse_data_new.py "$gff3_file_name" "$fasta_file_name" "${gff3_file_name%.gff3}.tsv"
rm "$gff3_file_name"
rm "$fasta_file_name"
rm *.gz