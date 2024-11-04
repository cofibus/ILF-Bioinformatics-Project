import requests
import json
from Bio import SeqIO
import pandas as pd
import pathlib


def oma_search(fasta_file):
    """Performs an OMA (Orthologous MAtrix) search for a given FASTA file.
    Results, including any errors, are stored in a JSON file.

    Args:
    fasta_file (str): The path to the input FASTA file.
    """
    with open(fasta_file, "r") as handle:
        sequences = []
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))

    base_url = "https://omabrowser.org/api/sequence/"
    params = {
        "search": "mixed",
        "full_length": "false"
    }

    all_results = {}  # Dictionary to store results for all sequences

    for seq in sequences:
        params["query"] = seq
        response = requests.get(base_url, params=params)

        if response.status_code == 200:
            all_results[seq] = response.json()  # Store the JSON response
        else:
            all_results[seq] = {
                "error": f"HTTP error {response.status_code}",
                "message": response.text
            }

    # Save the results to a JSON file
    with open("oma_search_results.json", "w") as output_file:
        json.dump(all_results, output_file, indent=4)  # Use indent for readability

def main():
  """Main function to execute the script."""
  full_dir = pathlib.Path(__file__).parent.absolute()
  fasta_file = full_dir / "ligand_sequences.fasta"
  oma_search(fasta_file)

if __name__ == "__main__":
  main()
