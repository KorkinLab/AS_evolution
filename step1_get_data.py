import pandas as pd
from ftplib import FTP
import os
import glob
import shutil

species_df = pd.DataFrame(columns=['database', 'mean_gene_length', 'variance_gene_length', 'mean_exon_count', 'variance_exon_count',  'mean_protein_length', 'variance_protein_length'])

if os.path.exists("temp_folder"):
    shutil.rmtree("temp_folder")
os.makedirs("temp_folder", exist_ok=True)
os.chdir("temp_folder")

ftp = FTP('ftp.ensemblgenomes.org')
ftp.login() 
ftp.cwd('/pub/release-61')

for kingdom in ['fungi', 'metazoa','plants', 'protists','plants']:
    print(f"Processing {kingdom}")
    os.makedirs(kingdom, exist_ok=True)

    ftp.cwd(f'/pub/release-61/{kingdom}/gff3')
    species_list = ftp.nlst()
    species_list = [s for s in species_list if 'collection' not in s]  # Filter out collection directories
    species_count = len(species_list)
    species_finished = 0
    printing_threshold = 0.1
    
    for species in species_list:
        # Remove all non-directory files in the current directory
        for f in os.listdir('.'):
            if os.path.isfile(f):
                os.remove(f)

        # Download the GFF3 file
        ftp.cwd(f'/pub/release-61/{kingdom}/gff3/{species}')
        gff_files = ftp.nlst()
        gff_files = sorted(gff_files, key=lambda f: ftp.size(f) if f.endswith('gff3.gz') else 0, reverse=True)
        gff_files = [f for f in gff_files if 'abinitio' not in f]
        filename = gff_files[0]
        if 'abinitio' in filename:
            print(f"Skipping abinitio file for {species}: {filename}")
            continue

        try:
            with open(filename, 'wb') as local_file_obj:
                ftp.retrbinary("RETR " + filename, local_file_obj.write)
        except Exception as e:
            print(f"Error downloading GFF3 for {filename}: {e}")
            continue


        # Download the corresponding FASTA file
        ftp.cwd(f'/pub/release-61/{kingdom}/fasta/{species}/pep')
        fa_files = ftp.nlst()
        fa_files = sorted(fa_files, key=lambda f: ftp.size(f) if f.endswith('fa.gz') else 0, reverse=True)
        fa_files = [f for f in fa_files if 'abinitio' not in f]
        filename = fa_files[0]

        try:
            with open(filename, 'wb') as local_file_obj:
                ftp.retrbinary("RETR " + filename, local_file_obj.write)
        except Exception as e:
            print(f"Error downloading fasta for {filename}: {e}")
            continue

        # Extract species name and database from the directory structure
        species_df.at[species, 'database'] = kingdom

        # Run the shell script to extract GTF data
        os.system(f"bash ../extract_data.sh")

        # Read the extracted data into the DataFrame
        tsv_files = glob.glob('*.tsv')
        if not tsv_files:
            print(f"No TSV files found for {species}. Skipping...")
            continue
        elif len(tsv_files) > 1:
            print(f"Multiple TSV files found for {species}. Using the first one: {tsv_files[0]}")
        tsv_file = tsv_files[0]
        extracted_data = pd.read_csv(tsv_file, sep='\t')

        # Calculate statistics
        species_df.at[species, 'mean_gene_length'] = extracted_data['gene_length'].mean()
        species_df.at[species, 'variance_gene_length'] = extracted_data['gene_length'].var()
        species_df.at[species, 'mean_exon_count'] = extracted_data['coding_exon_count'].mean()
        species_df.at[species, 'variance_exon_count'] = extracted_data['coding_exon_count'].var()
        species_df.at[species, 'mean_protein_length'] = extracted_data['protein_length'].mean()
        species_df.at[species, 'variance_protein_length'] = extracted_data['protein_length'].var()

        # Move the TSV file to the kingdom folder
        kingdom_folder = os.path.join(kingdom)
        new_tsv_path = os.path.join(kingdom_folder, os.path.basename(tsv_file))
        os.rename(tsv_file, new_tsv_path)

        species_finished += 1
        if species_finished >= species_count * printing_threshold:
            # Save species_df
            species_df.to_csv(os.path.join('../', 'panel_a_data.tsv'), sep='\t', index_label='gtf_name')
            print(f"Processed {printing_threshold*100:.1f}%")
            printing_threshold += 0.1


ftp.close()


ftp = FTP('ftp.ensembl.org')
ftp.login() 
ftp.cwd('/pub/release-114')

for kingdom in ['vertebrate']:
    print(f"Processing {kingdom}")
    os.makedirs(kingdom, exist_ok=True)
    
    ftp.cwd(f'/pub/release-114/gff3')
    species_list = ftp.nlst()
    species_count = len(species_list)
    species_finished = 0
    printing_threshold = 0.1
    
    for species in species_list:
        # Remove all non-directory files in the current directory
        for f in os.listdir('.'):
            if os.path.isfile(f):
                os.remove(f)

        # Download the GFF3 file
        ftp.cwd(f'/pub/release-114/gff3/{species}')
        gff_files = ftp.nlst()
        gff_files = sorted(gff_files, key=lambda f: ftp.size(f) if f.endswith('gff3.gz') else 0, reverse=True)
        gff_files = [f for f in gff_files if 'abinitio' not in f]
        filename = gff_files[0]
        if 'abinitio' in filename:
            print(f"Skipping abinitio file for {species}: {filename}")
            continue

        try:
            with open(filename, 'wb') as local_file_obj:
                ftp.retrbinary("RETR " + filename, local_file_obj.write)
        except Exception as e:
            print(f"Error downloading GFF3 for {filename}: {e}")
            continue

        # Download the corresponding FASTA file
        ftp.cwd(f'/pub/release-114/fasta/{species}/pep')
        fa_files = ftp.nlst()
        fa_files = sorted(fa_files, key=lambda f: ftp.size(f) if f.endswith('fa.gz') else 0, reverse=True)
        fa_files = [f for f in fa_files if 'abinitio' not in f]
        filename = fa_files[0]

        try:
            with open(filename, 'wb') as local_file_obj:
                ftp.retrbinary("RETR " + filename, local_file_obj.write)
        except Exception as e:
            print(f"Error downloading fasta for {filename}: {e}")
            continue

        # Extract species name and database from the directory structure
        species_df.at[species, 'database'] = kingdom

        # Run the shell script to extract GTF data
        os.system(f"bash ../extract_data.sh")

        # Read the extracted data into the DataFrame
        tsv_files = glob.glob('*.tsv')
        if not tsv_files:
            print(f"No TSV files found for {species}. Skipping...")
            continue
        elif len(tsv_files) > 1:
            print(f"Multiple TSV files found for {species}. Using the first one: {tsv_files[0]}")
        tsv_file = tsv_files[0]
        extracted_data = pd.read_csv(tsv_file, sep='\t')

        # Calculate statistics
        species_df.at[species, 'mean_gene_length'] = extracted_data['gene_length'].mean()
        species_df.at[species, 'variance_gene_length'] = extracted_data['gene_length'].var()
        species_df.at[species, 'mean_exon_count'] = extracted_data['coding_exon_count'].mean()
        species_df.at[species, 'variance_exon_count'] = extracted_data['coding_exon_count'].var()
        species_df.at[species, 'mean_protein_length'] = extracted_data['protein_length'].mean()
        species_df.at[species, 'variance_protein_length'] = extracted_data['protein_length'].var()

        # Move the TSV file to the kingdom folder
        kingdom_folder = os.path.join(kingdom)
        new_tsv_path = os.path.join(kingdom_folder, os.path.basename(tsv_file))
        os.rename(tsv_file, new_tsv_path)

        species_finished += 1
        species_df.to_csv(os.path.join('../', 'panel_a_data.tsv'), sep='\t', index_label='gtf_name')
        if species_finished >= species_count * printing_threshold:
            print(f"Processed {printing_threshold*100:.1f}%")
            printing_threshold += 0.1        
        

    ftp.cwd('..')  # Go back to the kingdom directory

ftp.close()
# Save the final DataFrame
species_df.to_csv(os.path.join('..', 'panel_a_data.tsv'), sep='\t', index_label='gtf_name')
# Clean up temporary folder
os.chdir('..')
os.rename("temp_folder", "ensembl_data")
