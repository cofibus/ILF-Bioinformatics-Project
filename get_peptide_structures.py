import pandas as pd 
import os
import requests
import pathlib
import Bio.PDB

def download_pdb_struct(pdb_id : str, save_path : pathlib.Path):
    """
    Downloads a PDB structure from the RCSB PDB database.
    """
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    r = requests.get(url)
    # check if it is not html
    if "the requested url was not found on this server" in r.text.lower():
        raise ValueError(f'{pdb_id} not found')

    with open(save_path, 'w') as f:
        f.write(r.text)


def get_chains_of_pdb(pdb_path : pathlib.Path):
    """
    Gets the chains of a PDB structure.
    """
    pdb_path = str(pdb_path)
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure('X', pdb_path)
    chains = structure[0].get_list()
    return chains


def keep_only_chain_of_pdb(pdb_path : pathlib.Path, chain_letter : str, save_path : pathlib.Path):
    """
    Keeps only a specific chain from a PDB structure and saves it to a new file.
    """
    pdb_path = str(pdb_path)
    save_path = str(save_path)

    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure('X', pdb_path)
    chain = structure[0][chain_letter]
    io = Bio.PDB.PDBIO()
    io.set_structure(chain)
    io.save(save_path)


def download_all():
    """
    Downloads and processes PDB structures.
    """
    script_dir = pathlib.Path(__file__).parent.absolute()
    full_dir = script_dir / 'hormone_complexes'
    full_dir.mkdir(exist_ok=True)
    out_dir = script_dir / 'hormone_structures_v2'
    out_dir.mkdir(exist_ok=True)

    # ,Peptide name,Sequence,Sequence length,PDB,Chain,Receptor,Family,GPCR Class
    df = pd.read_csv(script_dir / 'GPCRdb_peptide_ligands.csv')
    df['PDB'] = df['PDB'].str.upper()
    df['PDB'] = df['PDB'].apply(lambda x: x if x != '-' else '')
    df["Chain"] = df["Chain"].str.upper()
    df["Chain"] = df["Chain"].apply(lambda x: x if x != '-' else '')
    df["Chain"].fillna('', inplace=True)
    # fix filenames like '4XT1___CX<sub>3</sub>CL1.pdb''
    df["Peptide name"] = df["Peptide name"].str.replace('<sub>', '')
    df["Peptide name"] = df["Peptide name"].str.replace('</sub>', '')
    df["Peptide name"] = df["Peptide name"].str.replace('<sup>', '')
    df['Peptide name'] = df['Peptide name'].str.replace('</sup>', '')
    df["Peptide name"] = df["Peptide name"].str.replace('/', '_')

    for i, row in df.iterrows():
        # check if hormone has complex/structure
        pdb = row["PDB"]
        if len(pdb) == 0:
            continue
        # get structure info
        chain = row["Chain"]
        peptide = row["Peptide name"]
        # eliminating characters problematic for the foldseek search
        peptide = peptide.replace(" ", "_")
        peptide = peptide.replace(" ", "_")
        peptide = peptide.replace("&", "")
        peptide = peptide.replace("-", "_")
        peptide = peptide.replace(",", "_")
        peptide = peptide.replace(";", "_")
        peptide = peptide.replace("(", "_")
        peptide = peptide.replace(")", "_")
        peptide = peptide.replace("]", "_")
        peptide = peptide.replace("[", "_")

        # save complex
        complex_identifier = f'{pdb}___{peptide}'
        save_path_full = full_dir / f'{complex_identifier}.pdb'
        if save_path_full.exists():
            print(f'{save_path_full} exists')
        else:
            print(f'downloading {save_path_full}')
            try:
                download_pdb_struct(pdb, save_path_full)
            except ValueError as e:
                print(f'could not download {pdb}')
                print(e)
                continue
            print(f'downloaded {save_path_full}')

        # check chain
        chain_identifier = f'{peptide}_{pdb}'
        if len(chain) == 0:
            print(f'no chain for {save_path_full}')
            continue
        chains_in_pdb = get_chains_of_pdb(save_path_full)
        if chain.upper() not in [c.id.upper() for c in chains_in_pdb]:
            print(f'chain {chain} not in {save_path_full}')
            continue
        
        # save peptide structure
        save_path_chain = out_dir / f'{chain_identifier}.pdb'
        if save_path_chain.exists():
            print(f'{save_path_chain} exists')
        else:
            print(f'keeping chain {chain} of {save_path_full} to {save_path_chain}')
            keep_only_chain_of_pdb(save_path_full, chain, save_path_chain)
            print(f'kept chain {chain} of {save_path_full} to {save_path_chain}')


if __name__ == '__main__':
    download_all()
