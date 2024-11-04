#!/bin/bash -l
#SBATCH --qos=cpu
#SBATCH --partition=cpu_jobs                                         
#SBATCH --job-name foldseek_dwn
#SBATCH --mem=5G
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1             
#SBATCH --nodes=1   
#SBATCH --output=/projects/ilfgrid/data/lzr765/%j_out.txt
#SBATCH --error=/projects/ilfgrid/data/lzr765/%j_err.txt

echo "I am running a job"

BASE_DIR="/projects/ilfgrid/data/lzr765/foldseek_results"
if [[ ! -d "$BASE_DIR" ]]; then
    mkdir -p "$BASE_DIR"
fi

# Read the file containing IDs
while IFS=':' read -r name info; do
  # Extract the ID
  name="${name%%:*}"
  id=$(echo "$info" | jq -r '.id')

  if [[ -z "$id" ]]; then
    echo "Skipping line: Empty ID found"
    continue
  fi

  download_url="https://search.foldseek.com/api/result/download/$id"

  download_dir="$BASE_DIR/$name"
  if [[ ! -d "$download_dir" ]]; then
    mkdir -p "$download_dir"
  fi

  download_file="$download_dir/$name.tar.gz"

  curl -o "$download_file" "$download_url"

  if [[ $? -ne 0 ]]; then
    echo "Error downloading file for ID: $name"
    continue
  fi

  tar -xvzf "$download_file" -C "$download_dir"

  echo "Downloaded and unzipped file for ID: $name"
done < "/projects/ilfgrid/data/lzr765/foldseek_results.txt"

echo "I am done running a job"
# End of script