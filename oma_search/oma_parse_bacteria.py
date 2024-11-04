import json
import pathlib

def find_bacteria_targets(oma_data):
  """
  Extracts target information for bacterial species from the JSON data.

  Args:
    json_data: The JSON data as a Python dictionary.

  Returns:
    A list of dictionaries, where each dictionary contains information 
    about a bacterial target and the associated query.
  """
  bacteria_targets = []
  for query_key, query_data in oma_data.items():
    if 'targets' in query_data:
      for target in query_data['targets']:
        # Check if 'Bacteria' is in the hog_levels list
        if 'Bacteria' in target.get('hog_levels', []): 
          target_info = {
              'identified_by': query_data['identified_by'],
              'query': query_data['query'],
              'target_species': target['species']['species'],
              'target_sequence': target['sequence'],
              'target_description': target['description'],
              'hog_levels': target['hog_levels'], 
          }
          # Add alignment score and alignment only for approximate matches
          if query_data['identified_by'] == 'approximate match':
              target_info['alignment_score'] = target['alignment_score']
              target_info['alignment'] = target['alignment'][1]

          bacteria_targets.append(target_info)
  return bacteria_targets


def main():
  """Main function to execute the script."""
  full_dir = pathlib.Path(__file__).parent.absolute()
  json_file = full_dir / "oma_search_results.json"
  
  with open(json_file, 'r') as f:
    oma_data = json.load(f)
  
  bacteria_data = find_bacteria_targets(oma_data)
  print(len(bacteria_data))

  with open(full_dir / "oma_search_results_bacteria.json", "w") as output_file:
      json.dump(bacteria_data, output_file, indent=4)


if __name__ == "__main__":
  main()
