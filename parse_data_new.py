import sys
import pandas as pd

if len(sys.argv) != 4:
    print("Usage: python parse_gtf.py <gff3_file> <fasta_file> <tsv_file>")
    sys.exit(1)

gff3_file = sys.argv[1]
fasta_file = sys.argv[2]
tsv_file = sys.argv[3]

# parse gtf file
genes = []
gene_entry = [0,0,0,'']  # [gene_id, gene_length, coding_exon_count, Ensembl_canonical_transcript_id]
keep_reading = False

with open(gff3_file, "r") as file:
    for line in file:
        if line.startswith("#"):
            continue
        
        fields = line.strip().split("\t")

        feature = fields[2]

        attributes = fields[8].split(";")
        attributes_dict = {attr.split('=')[0].strip(): attr.split('=')[1].strip() for attr in attributes if '=' in attr}


        # check for protein coding genes
        if feature == 'gene':

            if attributes_dict.get('biotype', '') == "protein_coding":
                keep_reading = True
                genes.append(gene_entry)
                gene_entry = [0,0,0,'']  # [gene_id, gene_length, exon_count, Ensembl_canonical_transcript_id]
                canonical_transcript = ''
            else:
                keep_reading = False
                continue # skip non protein coding genes
        
        if keep_reading:
            start_pos = int(fields[3])
            stop_pos = int(fields[4])

            if feature == 'gene':                
                gene_entry[0] = attributes_dict.get('gene_id', '')  # gene_id
                gene_entry[1] = stop_pos - start_pos + 1  # gene length

            elif feature == "mRNA":                
                tags = attributes_dict.get('tag', '').split(',')
                if "Ensembl_canonical" in tags:
                    canonical_transcript = attributes_dict.get('transcript_id', '')
                    if 'version' in attributes_dict:
                        version = attributes_dict['version']
                        if canonical_transcript.split(".")[-1] != str(version):
                            gene_entry[3] = canonical_transcript + '.' + str(version)  # Ensembl_canonical_transcript_id
                        else:
                            gene_entry[3] = canonical_transcript
                    else:
                        gene_entry[3] = canonical_transcript  # Ensembl_canonical_transcript_id

            elif feature == "CDS":
                transcript_id = attributes_dict.get('Parent', '')[11:]               
                if transcript_id == canonical_transcript:
                    gene_entry[2] += 1



genes.append(gene_entry)  # append the last gene entry
# remove the first entry which is empty
genes = genes[1:]
# create a DataFrame and save to tsv
df = pd.DataFrame(genes, columns=["gene_id", "gene_length", "coding_exon_count", "Ensembl_canonical_transcript_id"])
df.set_index("Ensembl_canonical_transcript_id", inplace=True)
# df.to_csv(tsv_file, sep="\t", index=False)

canonical_transcript_list = df.index.tolist()
genes = []
gene_entry = ['',0]  # [transcript_id, protein_length]
# parse fasta file
keep_reading = True
with open(fasta_file, "r") as file:
    for line in file:
        if line.startswith(">"):
            transcript_id = line.strip(">").split(" ")[4][11:]  # extract transcript_id from the header
            if transcript_id not in canonical_transcript_list:
                keep_reading = False
                continue
            else:
                genes.append(gene_entry)
                gene_entry = [transcript_id, 0]  # reset gene_entry for the new transcript
                keep_reading = True
        else:
            if keep_reading:
                gene_entry[1] += len(line.strip())  # accumulate protein length

genes.append(gene_entry)  # append the last gene entry
# remove the first entry which is empty
genes = genes[1:]
df_protein = pd.DataFrame(genes, columns=["transcript_id", "protein_length"])
df_protein.set_index("transcript_id", inplace=True)
# merge the two DataFrames
concatenated_df = pd.concat([df, df_protein], axis=1)


if (concatenated_df['protein_length'] == 0).any():
    print(f"{tsv_file} contains 0 in protein_length")

elif concatenated_df['protein_length'].nunique() == 1:
    print(f"All values in 'protein_length' are the same: {tsv_file}")

concatenated_df.to_csv(tsv_file, sep="\t", index=False)