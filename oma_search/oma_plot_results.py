import json
import os
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import csv
from Bio import SeqIO
import random


def load_data(file_path):
  """
  Loads data from a JSON file.
  """
  with open(file_path, 'r') as f:
    return json.load(f)

def count_targets(bacteria_data):
  """
  Counts the number of targets for each query.
  """
  target_counts = {}
  for target_info in bacteria_data:
    query = target_info['query']
    target_counts[query] = target_counts.get(query, 0) + 1
  return target_counts

def get_query_names(queries_seq, fasta_file):
  """
  Gets the names of the queries from a FASTA file.
  """
  fasta_sequences = list(SeqIO.parse(fasta_file, "fasta"))
  query_names = []
  for query_seq in queries_seq:
    for record in fasta_sequences:
      if str(record.seq) == query_seq:
        # Keep both name and sequence
        record_name = record.description[:-2] 
        query_names.append(f"{record_name}\n{query_seq}") 
        break
    else:
      query_names.append(query_seq)  # If not found, keep the sequence
  return query_names

def create_bar_plot(bacteria_data, fasta_file, output_file):
  """
  Creates and saves a bar plot of target counts per query.
  """
  target_counts = count_targets(bacteria_data)

  queries = list(target_counts.keys())
  counts = list(target_counts.values())
  query_names = get_query_names(queries, fasta_file)

  plt.figure(figsize=(12, 8))
  
  bars = plt.bar(query_names, counts, color='orange')
  plt.xlabel("Query", fontsize='x-large')
  plt.ylabel("Number of bacterial targets", fontsize='x-large')
  plt.title("Bacterial target counts per query", fontsize='x-large')
  plt.xticks(rotation=45, ha="right", fontsize='large')

  for bar, count in zip(bars, counts):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), 
            str(count), ha='center', va='bottom')

  plt.tight_layout()

  plt.savefig(output_file)

def normalize_score(score, length):
  """
  Normalizes an alignment score by dividing it by the alignment length.
  """
  if length == 0:
    return 0.0
  return score / length


def create_scatter_plot(data, fasta_file, output_file):
  """
  Creates and saves a scatter plot of normalized alignment score vs. query length using Seaborn.
  """
  alignment_lengths = []
  normalized_scores = []
  queries_seq = []
  for target_info in data:
    aln_length = len(target_info['query']) if target_info['identified_by'] == 'exact match' else len(target_info['alignment'])
    alignment_lengths.append(aln_length)
    score = 100 if target_info['identified_by'] == 'exact match' else target_info['alignment_score']
    normalized_scores.append(normalize_score(score, aln_length))
    queries_seq.append(target_info['query'])

  plt.figure(figsize=(10, 6))
  sns.scatterplot(x=alignment_lengths, y=normalized_scores, hue=queries_seq, style=queries_seq, palette='tab20', markers=True, s=150)

  plt.xlabel("Alignment length", fontsize='x-large')
  plt.ylabel("Normalized alignment score", fontsize='x-large')
  plt.title("Normalized alignment score vs. alignment length", fontsize='x-large')

  unique_queries = list(set(queries_seq))
  query_names = get_query_names(unique_queries, fasta_file)

  query_to_name = dict(zip(unique_queries, query_names))

  handles, labels = plt.gca().get_legend_handles_labels()

  unique_labels = []
  unique_handles = []
  for label, handle in zip(labels, handles):
    if label not in unique_labels:
      unique_labels.append(label)
      unique_handles.append(handle)

  new_labels = [query_to_name[label].split("\n")[0] for label in unique_labels]  

  plt.legend(unique_handles, new_labels, title="Queries", loc='upper left', 
             bbox_to_anchor=(1, 1), fontsize='x-large')

  plt.tight_layout()
  plt.savefig(output_file, bbox_inches='tight')

def create_csv_report(bacteria_data, fasta_file, output_file):
  """
  Creates a CSV report of the bacterial target data.
  """
  fasta_sequences = list(SeqIO.parse(fasta_file, "fasta"))
  csv_data = []
  for target_info in bacteria_data:
    query_seq = target_info['query']
    for record in fasta_sequences:
      if str(record.seq) == query_seq:
        query_name = record.description[:-2]  # Remove the last two characters from the name
        break
    else:
      query_name = query_seq

    query_length = len(target_info['query'])
    target_length = len(target_info['target_sequence'])
    alignment = target_info.get('alignment')

    if target_info['identified_by'] == 'exact match':
      aligned_sequence = target_info['query']
      alignment_length = query_length
      alignment_start = target_info['target_sequence'].find(aligned_sequence) + 1 if aligned_sequence in target_info['target_sequence'] else -1
      alignment_end = alignment_start + alignment_length - 1 if alignment_start != -1 else -1
    else:
      aligned_sequence = alignment if alignment else "-"
      alignment_length = len(aligned_sequence)
      alignment_start = target_info['target_sequence'].find(aligned_sequence) + 1 if aligned_sequence in target_info['target_sequence'] else -1
      alignment_end = alignment_start + alignment_length - 1 if alignment_start != -1 else -1

    alignment_percentage = (alignment_length / target_length) * 100

    row = [
        target_info['target_species'],
        target_info['target_description'],
        query_name,
        target_info['identified_by'],
        aligned_sequence,
        target_length,
        f"{alignment_percentage:.2f}%",
        f"{alignment_start}-{alignment_end}"
    ]
    csv_data.append(row)

  headers = ["Target species", "Target description", "Query name", "Type of alignment",
             "Aligned sequence", "Target length", "Alignment percentage", "Aligned sequence position"]

  with open(output_file, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(headers)
    csv_writer.writerows(csv_data)

def main():
  """
  Main function to execute the script.
  """
  full_dir = pathlib.Path(__file__).parent.absolute()
  bacteria_file = full_dir / "oma_search_results_bacteria.json"
  fasta_seq_file = full_dir / "ligand_sequences.fasta"

  sns.set_style('whitegrid')
  sns.set_context('paper')

  plots_dir = full_dir / "oma_plots"
  if not plots_dir.exists():
    os.makedirs(plots_dir)

  bacteria_data = load_data(bacteria_file)

  create_bar_plot(bacteria_data, fasta_seq_file, plots_dir / "bacteria_target_counts.png")
  create_scatter_plot(bacteria_data, fasta_seq_file, plots_dir / "alignment_score_vs_length.png")
  create_csv_report(bacteria_data, fasta_seq_file, plots_dir / "bacteria_targets_table.csv")


if __name__ == "__main__":
  main()
