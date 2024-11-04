import pandas as pd
import numpy as np
import os
import sys
import pathlib
import time

from species_name_to_taxon_id import get_taxon_id_uniprot
from taxon_to_lineage import get_taxon_lineage_batch


def get_df():
    """ 
    Read dataframe

    Header: 
    Query,Filename,Job ID,Target/Description,pident,alnlen,mismatch(not sure),gapopen,qstart,Qend,tstart,tend,Prob.,
    E-value,Score(not sure),Query length,Target length,Qaln,Taln,tca,tseq,Taxid,taxname/species
    """
    script_dir = pathlib.Path(__file__).parent.absolute()
    df_p = script_dir / 'foldseek_parsed_results_nohuman.csv'
    df = pd.read_csv(df_p)
    return df


def get_species_to_taxon(species, species_csv, sleep_time=0.1, save_interval=100):
    """
    Get taxons for species in the dataframe and save to a csv file
    """
    # load or create the species to taxon csv
    if os.path.exists(species_csv):
        species_to_taxon = pd.read_csv(species_csv, index_col=0)
    else:
        species_to_taxon = pd.DataFrame(columns=['Species', 'Taxon ID'])
        species_to_taxon = species_to_taxon.set_index('Species')

    species_found = set(species_to_taxon.index)
    species_missing = set(species) - species_found

    # expected runtime
    print(f"Expected runtime: {round(sleep_time * len(species_missing))} seconds")
    for i, species in enumerate(species_missing):
        taxon_id = get_taxon_id_uniprot(species)
        species_to_taxon.loc[species] = taxon_id
        if i % save_interval == 0:
            species_to_taxon.to_csv(species_csv)
        time.sleep(sleep_time)

    species_to_taxon['Taxon ID'] = species_to_taxon['Taxon ID'].apply(lambda x: np.nan if 'Error' in str(x) else x)    
    species_to_taxon['Taxon ID'] = species_to_taxon['Taxon ID'].astype('Int64')
    species_to_taxon.to_csv(species_csv)
    return species_to_taxon


def get_taxon_to_lineage(taxon_ids, taxon_csv, email,
                         batch_size=50, sleep_time=0.1, save_interval=100):
    """ 
    Get taxon to lineage and save to a csv file
    """
    # load or create the taxon to lineage csv
    if os.path.exists(taxon_csv):
        taxon_to_lineage = pd.read_csv(taxon_csv, index_col=0)
    else:
        taxon_to_lineage = pd.DataFrame(columns=['Taxon ID', 'Lineage'])
        taxon_to_lineage = taxon_to_lineage.set_index('Taxon ID')

    taxon_ids_found = set(taxon_to_lineage.index)
    taxon_ids_missing = set(taxon_ids) - taxon_ids_found
    
    # split missing into batches
    if len(taxon_ids_missing) > 0:
        batch_count = len(taxon_ids_missing) // batch_size + 1
        print("Expected runtime: ", round(sleep_time * len(taxon_ids_missing) * batch_count), "seconds")
    else:
        batch_count = 0

    batches = []
    for i in range(batch_count):
        batch = list(taxon_ids_missing)[i*batch_size:(i+1)*batch_size]
        batches.append(batch)
   
    for i, batch in enumerate(batches):
        lineages = get_taxon_lineage_batch(batch, email)
        for taxon_id, lineage in lineages.items():
            taxon_to_lineage.loc[taxon_id] = lineage
        if i % save_interval == 0:
            taxon_to_lineage.to_csv(taxon_csv)
        time.sleep(sleep_time)
    
    # make sure taxon id is int
    taxon_to_lineage.index = taxon_to_lineage.index.astype(int)

    taxon_to_lineage.to_csv(taxon_csv)
    return taxon_to_lineage


def get_phylogentic_options():
    """
    List of phylogenetic classification levels.
    """
    return ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def parse_lineage_string(l:str) -> dict:
    """
    Parse lineage string into a dictionary of phylogenetic levels.
    """
    options = get_phylogentic_options()
    parsed = {}
    split_parts = l.split(';')
    split_parts = [part.strip() for part in split_parts]

    for option_i, option in enumerate(options):
        if option_i >= len(split_parts):
            parsed[option] = np.nan
        else:
            parsed[option] = split_parts[option_i]
    return parsed


def get_lineage_part(lineage :str, part :str) -> str:
    """
    Extract a specific part (e.g., Domain, Phylum) from a lineage string.

    d__ for Domain
    p__ for Phylum
    c__ for Class
    o__ for Order
    f__ for Family
    g__ for Genus
    s__ for Species
    """
    options = get_phylogentic_options()
    assert part in options, f"part must be one of {options}"
    if pd.isna(lineage):
        return None
    parsed = parse_lineage_string(lineage)
    return parsed[part]


def get_lineage(df, id, ignore_nan=True):
    """
    Retrieve the lineage string for a given ID from a DataFrame.
    """
    if id not in df.index:
        return None
    else:
        lineage = df.loc[id, 'Lineage']
        if ignore_nan and pd.isna(lineage):
            return None
        return lineage


def main():
    df = get_df()
    script_dir = pathlib.Path(__file__).parent.absolute()
    plot_dir = script_dir / 'plots'
    plot_dir.mkdir(parents=True, exist_ok=True)

    print("Number of rows: ", len(df))
    print(df['Filename'].value_counts())
    print("Number of unique species", len(df['taxname/species'].unique()))

    # get species to taxon
    species = df['taxname/species'].unique()
    species_to_taxon = get_species_to_taxon(species, script_dir / 'foldseek_species_to_taxon.csv',
                                            sleep_time=2, save_interval=50)
    
    # add taxon id to df
    df['Taxon ID'] = df['taxname/species'].apply(lambda x: species_to_taxon.loc[x])
    taxon_ids = df['Taxon ID'].dropna().unique()

    # save taxons    
    print("Number of unique taxons", len(df['Taxon ID'].unique()))
    taxon_p = script_dir / 'foldseek_taxids.txt'
    taxon_p.write_text('\n'.join([str(int(t)) for t in taxon_ids]))
    
    # get taxon to lineage
    taxon_to_lineage = get_taxon_to_lineage(taxon_ids, script_dir / 'foldseek_taxon_to_lineage.csv',
                                            'lzr765@ku.dk', batch_size=100,
                                            sleep_time=1, save_interval=5)
    
    # add lineage to df (ignore nans)
    df['Lineage'] = df['Taxon ID'].apply(lambda x: get_lineage(taxon_to_lineage, x, ignore_nan=True))

    # add each individual part of the lineage to the df
    options = get_phylogentic_options()
    for option in options:
        df[option] = df['Lineage'].apply(get_lineage_part, part=option)

    # def check if is in taxon list from GMrepo
    GMrepo_taxon_p = script_dir / 'GMrepo_species_taxon_ids_morethan3.txt'
    known_gut_microbe_taxons = set(map(int, GMrepo_taxon_p.read_text().split('\n')))
    df['KnownGutMicrobe(GMrepo)'] = df['Taxon ID'].apply(lambda x: x in known_gut_microbe_taxons)

    print(df['KnownGutMicrobe(GMrepo)'].value_counts())
    print(df['kingdom'].value_counts())

    # save the df
    df.to_csv(script_dir / 'foldseek_parsed_results_nohuman_annotated.csv')

if __name__ == "__main__":
    main()
