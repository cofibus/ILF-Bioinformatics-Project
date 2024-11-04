import requests

def get_taxon_id_uniprot(species_name : str) -> int:
    """
    Retrieve the taxon ID for a given species name from UniProt.
    """
    url = "https://rest.uniprot.org/taxonomy/search"
    params = {
        "query": species_name,
        "format": "json"
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.json()
        if "results" in data and len(data["results"]) > 0:
            taxon_id = data["results"][0]["taxonId"]
            return taxon_id
        else:
            return None
    else:
        return f"Error: {response.status_code} - {response.text}"


def check_get_taxon_id_uniprot():
    """
    Test the get_taxon_id_uniprot function with known species names.
    """
    assert get_taxon_id_uniprot('Atractaspis engaddensis') == 1343144, 'Atractaspis engaddensis error'
    assert get_taxon_id_uniprot('synthetic construct') == 32630, 'synthetic construct error'
    assert get_taxon_id_uniprot('Platichthys flesus') == 8260, 'Platichthys flesus error'
    print("All tests passed for check_get_taxon_id_uniprot")

if __name__ == "__main__":
    # Example
    print(get_taxon_id_uniprot('Atractaspis engaddensis'))
