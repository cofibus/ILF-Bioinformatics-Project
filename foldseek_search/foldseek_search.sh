#!/bin/bash -l
#SBATCH --qos=cpu
#SBATCH --partition=cpu_jobs                                         
#SBATCH --job-name foldseek_search
#SBATCH --mem=5G
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1             
#SBATCH --nodes=1   
#SBATCH --output=/projects/ilfgrid/data/lzr765/%j_out.txt
#SBATCH --error=/projects/ilfgrid/data/lzr765/%j_err.txt

# Performs Foldseek searches for all PDB files in a directory.

# It iterates through each PDB file, extracts the peptide name from the filename,
# constructs the Foldseek API command using curl, and appends the peptide name and
# the search results (ticket ID) to an output file.

# It also includes a delay between each search to avoid rate limiting.

pdb_dir="/projects/ilfgrid/data/lzr765/computer_files/peptide_structures"

output_file="/projects/ilfgrid/data/lzr765/foldseek_results.txt"

for pdb_file in "$pdb_dir"/*.pdb; do
  echo "Processing: $pdb_file"

  peptide_name=$(basename "$pdb_file" .pdb)

  curl_output=$(curl -X POST \
    -F "q=@$pdb_file" \
    -F 'mode=3diaa' \
    -F 'database[]=afdb50' \
    -F 'database[]=afdb-swissprot' \
    -F 'database[]=afdb-proteome' \
    -F 'database[]=cath50' \
    -F 'database[]=mgnify_esm30' \
    -F 'database[]=pdb100' \
    -F 'database[]=gmgcl_id' \
    -F 'database[]=bfmd' \
    https://search.foldseek.com/api/ticket)


  echo "$peptide_name: $curl_output" >> "$output_file"

  echo "Finished: $pdb_file"

  sleep 180
done
# End of script
