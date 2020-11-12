"""get_uniprot_features.py obtains info from uniprot for the human genome and
individual uniprot ids

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-08-21

"""

# Standard library imports
import urllib.request

# Third party imports
import pandas as pd


def extract_uniprot_info(uniprot_id, reviewed='*'):
    """Uniprot search request
    Get info from here (https://www.uniprot.org/help/uniprotkb_column_names)
    """
    search_string = (
        'https://www.uniprot.org/uniprot/'
        f'?query=reviewed:{reviewed}'
        f'+AND+accession:{uniprot_id}'
        '+&format=tab&columns=id,entry%20name,protein%20name,'
        'genes(PREFERRED),organism,sequence,database(PDB),database(DisProt)'
        ',database(MobiDB),length,feature(TRANSMEMBRANE),'
        'feature(INTRAMEMBRANE),feature(TOPOLOGICAL%20DOMAIN),'
        'comment(SIMILARITY),'
        'comment(DOMAIN),'
        'feature(DOMAIN%20EXTENT),'
        'families,'
        'feature(COILED%20COIL),'
        'feature(MOTIF),'
        'feature(REGION),'
        'feature(REPEAT),'
        'feature(ZINC%20FINGER),'
        'feature(DISULFIDE%20BOND),'
        'feature(ACTIVE%20SITE),'
        'feature(BINDING%20SITE),'
        'database(Pfam),'
        'database(InterPro),'
        'database(SUPFAM),'
        'database(PROSITE),'
        'database(MUTAGENESIS),'
        'comment(ALLERGEN),'
        'comment(BIOTECHNOLOGY),'
        'comment(DISRUPTION%20PHENOTYPE),'
        'comment(DISEASE),'
        'comment(PHARMACEUTICAL),'
        'comment(TOXIC%20DOSE),'
        'feature(NATURAL%20VARIANT),'
        'feature(GLYCOSYLATION),'
        'feature(INITIATOR%20METHIONINE),'
        'feature(LIPIDATION),'
        'feature(MODIFIED%20RESIDUE),'
    )
    req2 = urllib.request.Request(search_string)
    with urllib.request.urlopen(req2) as f:
       response2 = f.read()
    result = [i.split("\t") for i in response2.decode("utf-8").split("\n")]
    result = pd.DataFrame(data=result[1:], columns=result[0])
    
    return result


def extract_uniprot_human_proteome_info(reviewed='*'):
    """Uniprot search request for human proteome"""
    search_string = (
        'https://www.uniprot.org/uniprot/'
        f'?query=reviewed:{reviewed}'
        '+AND+proteome:up000005640'
        '+&format=tab&columns=id,entry%20name,protein%20name,'
        'genes(PREFERRED),organism,sequence,database(PDB),database(DisProt)'
        ',database(MobiDB),length,feature(TRANSMEMBRANE),'
        'feature(INTRAMEMBRANE),feature(TOPOLOGICAL%20DOMAIN),'
        'comment(SIMILARITY),'
        'comment(DOMAIN),'
        'feature(DOMAIN%20EXTENT),'
        'families,'
        'feature(COILED%20COIL),'
        'feature(MOTIF),'
        'feature(REGION),'
        'feature(REPEAT),'
        'feature(ZINC%20FINGER),'
        'feature(DISULFIDE%20BOND),'
        'feature(ACTIVE%20SITE),'
        'feature(BINDING%20SITE),'
        'database(Pfam),'
        'database(InterPro),'
        'database(SUPFAM),'
        'database(PROSITE),'
    )
    req2 = urllib.request.Request(search_string)
    with urllib.request.urlopen(req2) as f:
       response2 = f.read()
    result = [i.split("\t") for i in response2.decode("utf-8").split("\n")]
    result = pd.DataFrame(data=result[1:], columns=result[0])
    
    return result


def extract_by_uniprot_fasta(keyword):
    """extract information from uniprot"""
    url_base = "https://www.uniprot.org/uniprot/"
    search_params = "?query=reviewed:yes" +\
        "+AND+accession:" + keyword
    return_params = "+&format=tab&columns=id,sequence,organism,entry%20name"
    url = url_base + search_params + return_params

    data_array = []
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8')
        unprocessed = line[:-1].split('\t')
        data_array.append(unprocessed)

    return data_array


def uniprot_search(key, args):
    """
    Uniprot search request, adapted from XXXXX
    """
    import urllib.parse
    import urllib.request

    if args.gene:
        keys = "+AND+gene:" + key
    else:
        keys = "+AND+"+key

    if args.organism:
        organism = "+AND+organism:" + args.organism
    else:
        organism = ""

    req2 = urllib.request.Request("https://www.uniprot.org/uniprot/?query=reviewed:yes"+organism+keys+"+&format=tab&columns=id,entry%20name,protein%20name,genes(PREFERRED),organism,database(PDB),database(DisProt),database(MobiDB),length,feature(TRANSMEMBRANE),feature(INTRAMEMBRANE),feature(TOPOLOGICAL%20DOMAIN),comment(SIMILARITY)"+args.additional)
    if args.uniprot:
        reg2 = urllib.request.Request("https://www.uniprot.org/uniprot/?query=reviewed:yes"+keys+"+&format=tab&columns=id,entry%20name,protein%20name,genes(PREFERRED),organism,database(PDB),database(DisProt),database(MobiDB),length,feature(TRANSMEMBRANE),feature(INTRAMEMBRANE),feature(TOPOLOGICAL%20DOMAIN),comment(SIMILARITY)"+args.additional)
    with urllib.request.urlopen(req2) as f:
       response2 = f.read()
    result = [i.split("\t") for i in response2.decode("utf-8").split("\n")]

    return result
