""" Show plots for the foldseek results for human hormones.
"""
import pathlib
import os 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd 


def get_df():
    """
    Load the annotated foldseek results.

    Header:
    Row,Query,Filename,Job ID,Target/Description,pident,alnlen,mismatch(not sure),gapopen,qstart,Qend,tstart,tend,Prob.,
    E-value,Score(not sure),Not sure,Target length,Qaln,Taln,tca,tseq,Taxid,taxname/species
    """
    #     
    script_dir = pathlib.Path(__file__).parent.absolute()
    df_p = script_dir / 'foldseek_parsed_results_nohuman_annotated.csv'
    df = pd.read_csv(df_p)

    df = df.rename(columns={
        # ... other column mappings
        'KnownGutMicrobe(GMrepo)': 'Known Gut Microbe',
        'Filename': 'Database'
    })

    filename_mapping = {
    'alis_pdb100.m8': 'PDB100',
    'alis_afdb-swissprot.m8': 'AFDB-Swissprot',
    'alis_afdb50.m8': 'AFDB50',
    'alis_mgnify_esm30.m8': 'MGnify_ESM30',
    'alis_afdb-proteome.m8': 'AFDB-Proteome',
    'alis_gmgcl_id.m8': 'GMGCL_ID',
    'alis_bfmd.m8': 'BFMD',
    'alis_cath50.m8': 'CATH50'
    }
    
    df['Database'] = df['Database'].map(filename_mapping)

    return df

def save(save_p):
    """
    Save plots in SVG and PNG formats.
    """
    plt.savefig(save_p)
    plt.savefig(save_p.with_suffix('.png'), dpi=300)
    plt.close()


def plot_kingdoms(df, save_p):
    """
    Plot the distribution of kingdoms.
    """
    kingdom_counts = df['Kingdom'].value_counts()
    plt.figure(figsize=(5, 5))
    ax = sns.barplot(x=kingdom_counts.values, y=kingdom_counts.index, log=True)

    for p in ax.patches:
        width = int(p.get_width())
        scientific_label = f'{width:,}'
        ax.annotate(scientific_label,
                    (p.get_width(), p.get_y() + p.get_height() / 2),
                    ha='left', va='center')

    plt.ylabel('')
    plt.xlabel('Count')
    plt.xlim(None, kingdom_counts.values.max() * 10)
    plt.title('Kingdom distribution of Foldseek\nresults (all databases)') 
    plt.tight_layout()
    save(save_p)


def plot_databases(df, save_p):
    """
    Plot which database of origin the hits are from.
    Only bacteria
    """
    df = df[df['Kingdom'] == 'Bacteria']
    # Make filename categorical
    df['Database'] = pd.Categorical(df['Database'])

    plt.figure(figsize=(5, 5))
    ax = sns.countplot(y='Database', data=df, hue='Known Gut Microbe', palette='magma',
                          log=True, order=df['Database'].value_counts().index)
    
    for p in ax.patches:
        width = int(p.get_width())
        scientific_label = f'{width:,}'
        ax.annotate(scientific_label,
                    (p.get_width(), p.get_y() + p.get_height() / 2),
                    ha='left', va='center')
        
    plt.ylabel('')
    plt.xlabel('Count')
    plt.xlim(None, df['Database'].value_counts().values.max() * 10)
    plt.title('Database distribution of Foldseek\nresults (Bacteria)')
    plt.tight_layout()
    save(save_p)


def plot_query_len_target_len(df, save_p, hue='pident', hue_threshold=None):
    """ 
    Plot identity score vs target length.
    Only known bacteria from gut microbiome
    """
    df = df[df['Kingdom'] == 'Bacteria']
    df = df[df['Known Gut Microbe']]
    df['pident'] = df['pident'].apply(float)
    
    # if hue threshold, remove anything BELOW threshold
    if hue_threshold:
        df = df[df[hue] >= hue_threshold]
        
    plt.figure(figsize=(8, 5))
    
    cmap = sns.color_palette("rainbow", as_cmap=True)
    df_filt = df[df['Query length'] <= 150]
    top_5 = df_filt.nlargest(5, hue)

    ax = sns.scatterplot(x='Target length', y='Query length', data=df,
                         alpha=0.5, hue=hue, palette=cmap, legend=False, edgecolor='black')

    #ax = sns.scatterplot(x='Target length', y='Query length', data=top_5, alpha=1, color='orange', edgecolor='black', s=40, ax=ax)
    
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                 norm=plt.Normalize(vmin=df[hue].min(), vmax=df[hue].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, shrink=0.3, ax=ax)         
    cbar.set_label(hue)

    plt.xlim(0, None)
    plt.ylim(0, 500)

    plt.xlabel('Query length')
    plt.ylabel('Target length')
    plt.title('Query length vs Target length\nof known gut microbes')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout(pad=2.0)

    save(save_p)
    return top_5

def plot_query_len_target_len_evalue(df, save_p, hue='E-value', hue_threshold=None):
    """ 
    Plot identity score vs target length.
    Only known bacteria from gut microbiome
    """
    df = df[df['Kingdom'] == 'Bacteria']
    df = df[df['Known Gut Microbe']]
    df['pident'] = df['pident'].apply(float)
    
    # if hue threshold, remove anything BELOW threshold
    if hue_threshold:
        df = df[df[hue] >= hue_threshold]
        
    plt.figure(figsize=(8, 5))
    
    cmap = sns.color_palette("rainbow", as_cmap=True)
    df_filt = df[df['Query length'] <= 150]
    top_5 = df_filt.nsmallest(5, hue)

    cmap = cmap.reversed()  # so that colors match (good evalues and good identity scores)

    ax = sns.scatterplot(x='Target length', y='Query length', data=df,
                         alpha=0.5, hue=hue, palette=cmap, legend=False, edgecolor='black')

    #ax = sns.scatterplot(x='Target length', y='Query length', data=top_5, alpha=1, color='orange', edgecolor='black', s=40, ax=ax)
    
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                 norm=plt.Normalize(vmin=df[hue].min(), vmax=df[hue].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, shrink=0.3, ax=ax)         
    cbar.set_label(hue)

    plt.xlim(0, None)
    plt.ylim(0, 500)

    plt.xlabel('Query length')
    plt.ylabel('Target length')
    plt.title('Query length vs Target length\nof known gut microbes')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()

    save(save_p)
    return top_5

def find_strong_matches(df):
    df = df[df['Kingdom'] == 'Bacteria']
    df = df[df['Known Gut Microbe']]
    #df_filt = df[df['Query length'] <= 150]

    strong_matches_df = df[(df['pident'].astype(float) >= 10) & (df['Query length'] <= 400)]
    peptide_counts = strong_matches_df['Query'].value_counts()
    #peptide_counts = df_filt['Query'].value_counts()
    print("Peptide strong matches counts:")
    print(peptide_counts.head(5))

    # Sort by pident in descending order and then by E-value in ascending order
    top5_both = strong_matches_df.sort_values(['pident', 'E-value'], ascending=[False, True]).head(5)

    return top5_both

def main():
    sns.set_style('whitegrid')
    sns.set_context('paper')

    df = get_df()
    #df_filt = df[(df['Target length'] <= 500) & (df['Query length'] <= 150)]

    script_dir = pathlib.Path(__file__).parent.absolute()
    plot_dir = script_dir / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    #databases = pd.unique(df['Database'])
    #print(databases)

    plot_kingdoms(df.copy(), plot_dir/ 'kingdoms.svg')
    plot_databases(df.copy(), plot_dir / 'databases.svg')
    
    top_both = find_strong_matches(df.copy())
    plot_query_len_target_len(df.copy(), plot_dir / 'query_target_identity.svg', hue='pident')
    top_identity20 = plot_query_len_target_len(df.copy(), plot_dir / 'query_target_identity20.svg', hue='pident', hue_threshold=20)
    top_evalue = plot_query_len_target_len_evalue(df.copy(), plot_dir / 'query_target_evalue.svg', hue='E-value') 
    top_prob = plot_query_len_target_len(df.copy(), plot_dir / 'query_target_prob.svg', hue='Prob.')

    top_data = [
        ("top_identity20", top_identity20),
        ("top_evalue", top_evalue),
        ("top_prob", top_prob),
        ("top_both", top_both)
    ]

    csv_dir = script_dir / 'top_csv'
    csv_dir.mkdir(exist_ok=True)
    combined_df = pd.concat([df.assign(key=name) for name, df in top_data])

    # Save the combined DataFrame to a single CSV file
    csv_path = csv_dir / "top_all.csv"
    combined_df.to_csv(csv_path, index=False)
        


if __name__ == "__main__":
    main()
