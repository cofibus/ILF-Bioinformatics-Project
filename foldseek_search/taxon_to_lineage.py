from Bio import Entrez

def get_taxon_lineage_batch(taxon_ids : list, email : str):
    """
    Retrieve the taxonomic lineage for a batch of taxon IDs using Entrez.
    """
    Entrez.email = email
    lineages = {}

    # Convert the list of IDs into a comma-separated string
    ids_str = ','.join(map(str, taxon_ids))

    try:
        # Fetch taxonomy data for all IDs in one request
        handle = Entrez.efetch(db="taxonomy", id=ids_str, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        # Extract lineage for each ID
        for r_i, record in enumerate(records):
            # taxon_id = record['TaxId']
            taxon_id = taxon_ids[r_i]
            lineage = record['Lineage']
            lineages[taxon_id] = lineage

    # try one by one, and return nan if keeps failing
    except Exception as e:
        print(f"Error retrieving Taxon IDs: {e}")
        for taxon_id in taxon_ids:
            try:
                handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
                record = Entrez.read(handle)
                handle.close()
                lineage = record[0]['Lineage']
                lineages[taxon_id] = lineage
            except Exception as e:
                print(f"Error retrieving Taxon ID {taxon_id}: {e}")
                lineages[taxon_id] = None

    return lineages

if __name__ == "__main__":
    # Example batch
    lineages = get_taxon_lineage_batch(["9606", "544730", "708627"]*50, 'lzr765@ku.dk')
    print(lineages)
