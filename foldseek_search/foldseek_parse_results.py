import os
import pandas as pd

def parse_foldseek_results(results_dir):
    """
    Parses FoldSeek results files from each folder into a DataFrame. Filtered by non-human results.
    """

    data = []
    for subdir, _, files in os.walk(results_dir):
        for file in files:
            if file.endswith('.m8'):
                file_path = os.path.join(subdir, file)
                with open(file_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line.startswith('#'):
                            columns = line.split('\t')
                            if len(columns) >= 6:
                                directory_name = os.path.basename(subdir) 
                                data.append([directory_name] + [file] + columns)


    # Column names based on the provided example
    column_names = [
        'Query', 'Filename', 'Job ID', 'Target and Description', 'pident', 'alnlen', 'Mismatch(not sure)',
        'Gapopen', 'Query start', 'Query end', 'Target start', 'Target end',
        'Probability', 'E-value', 'Score(not sure)', 'Query length',
        'Target length', 'Qaln', 'Taln', 'tca', 'tseq', 'Taxid', 'taxname/species'
    ]
    #Query,Filename,Job ID,Target/Description,pident,alnlen,mismatch(not sure),gapopen,qstart,Qend,tstart,tend,Prob.,
    # E-value,Score(not sure),Not sure,Target length,Qaln,Taln,tca,tseq,Taxid,taxname/species
   

    df = pd.DataFrame(data, columns=column_names)
    df = df[df.iloc[:, -1] != "Homo sapiens"]
    return df

def main():
  """
  Main function to execute the script.
  """
  full_dir = "/projects/ilfgrid/data/lzr765"
  results_dir = full_dir / "foldseek_results"
  df = parse_foldseek_results(results_dir)
  df.to_csv(full_dir / "foldseek_parsed_results_nohuman.csv", index=False)

if __name__ == "__main__":
  main()
