import pandas as pd
import json
import pathlib

def load_data(full_dir):
  """Loads peptide data from a JSON file and an Excel file.

  Args:
    full_dir Path to the directory containing the necessary files.
    
  Returns:
    A tuple containing:
      - A list of peptide dictionaries.
      - A list of unique ligands.
  """
  info_file = full_dir / "GPCRdb_peptide_ligands_info.txt"
  peptide_list = full_dir / "GPCRdb_peptides.xls"

  with open(info_file, 'r') as file:
      peptide_dicts = json.load(file)
  
  df = pd.read_excel(peptide_list)

  unique_ligands = df['ligand'].dropna().unique()
  unique_ligands = [ligand.strip().lower() for ligand in unique_ligands]
  return peptide_dicts, unique_ligands

def create_ligand_sequence_dict(peptide_dicts):
  """Creates a dictionary mapping ligand names to lists of sequences.

  Args:
    peptide_dicts: A list of peptide dictionaries.

  Returns:
    A dictionary where keys are ligand names and values are lists of sequences.
  """
  all_sequences_dict = {}
  for d in peptide_dicts:
      ligand = d['Peptide name'].strip().lower()
      sequence = d['Sequence']
      if sequence is not None:
          if ligand in all_sequences_dict:
              all_sequences_dict[ligand].append(sequence)
          else:
              all_sequences_dict[ligand] = [sequence]
  return all_sequences_dict

def create_correction_dict():
  """Creates a dictionary for correcting ligand names.

  Returns:
    A dictionary mapping incorrect ligand names to correct ones.
  """
  correct_id = {
      "cxcl12α": "cxcl12&alpha;",
      "[pyr1]apelin-13": "[pyr<sup>1</sup>]apelin-13",
      "[des-arg10]kallidin": "[des-arg<sup>10</sup>]kallidin",
      "[des-arg9]bradykinin": "[des-arg<sup>9</sup>]bradykinin",
      "[hyp3]bradykinin": "[hyp<sup>3</sup>]bradykinin",
      "lys-[hyp3]-bradykinin": "lys-[hyp<sup>3</sup>]-bradykinin", 
      "α-cgrp": "&alpha;-cgrp",  
      "β-cgrp": "&beta;-cgrp", 
      "cxcl12γ": "cxcl12&gamma;",  
      "cxcl12δ": "cxcl12&delta;", 
      "cxcl12φ": "cxcl12&phi;", 
      "cxcl12ε": "cxcl12&epsilon;", 
      "cxcl12β": "cxcl12&beta;", 
      "prp106-126": "prp<sub>106-126</sub>",
      "[des-gln14]ghrelin": "[des-gln<sup>14</sup>]ghrelin", 
      "α-msh": "&alpha;-msh", 
      "β-msh": "&beta;-msh", 
      "γ-msh": "&gamma;-msh", 
      "β-endorphin": "&beta;-endorphin", 
      "α-neoendorphin": "&alpha;-neoendorphin", 
      "β-neoendorphin": "&beta;-neoendorphin", 
      "neuropeptide-γ": "neuropeptide-&gamma;", 
      "neuropeptide γ": "neuropeptide &gamma;", 
      "prokineticin-2β": "prokineticin-2&beta;"
  }
  return correct_id

def find_sequences(peptide_dicts, unique_ligands):
  """Finds matching sequences for ligands.

  Args:
    peptide_dicts: A dictionary mapping ligand names to lists of sequences.
    unique_ligands: A list of peptide ligands.

  Returns:
    A tuple containing:
      - A dictionary mapping ligands to their sequences.
      - A list of ligands with missing sequences.
  """

  all_sequences_dict = create_ligand_sequence_dict(peptide_dicts)

  correct_id = create_correction_dict()

  final_peptides_and_seq = {}
  missing_seq = []

  for ligand in unique_ligands:
      seq = all_sequences_dict.get(ligand, None)
      if seq is None and ligand in correct_id:
          corrected_ligand = correct_id[ligand]
          seq = all_sequences_dict.get(corrected_ligand, None)
      if seq is not None:
          final_peptides_and_seq[ligand] = seq
      else:
          missing_seq.append(ligand)
  return final_peptides_and_seq, missing_seq

def write_fasta(final_peptides_and_seq, full_dir):
  """Writes the results to a FASTA file.

  Args:
    final_peptides_and_seq: A dictionary mapping ligands to their sequences.
    filename: The path to the output FASTA file.
  """

  filename = full_dir / "ligand_sequences.fasta"
  with open(filename, "w") as fasta_file:
      for ligand, sequences in final_peptides_and_seq.items():
          for idx, sequence in enumerate(sequences):
              fasta_file.write(f">{ligand}_{idx + 1}\n")
              fasta_file.write(f"{sequence}\n")

def main():
  """Main function to execute the script."""
  full_dir = pathlib.Path(__file__).parent.absolute()
  peptide_dicts, unique_ligands = load_data(full_dir)
  print(len(unique_ligands))
  final_peptides_and_seq, missing_seq = find_sequences(peptide_dicts, unique_ligands)
  write_fasta(final_peptides_and_seq, full_dir)

if __name__ == "__main__":
  main()
