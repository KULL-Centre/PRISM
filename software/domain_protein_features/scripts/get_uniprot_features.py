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
    
    features = [
		#standard features:
		'genes(PREFERRED)','organism','length','sequence',
		#database crossrefs
		'database(PDB)','database(DisProt)',
		'database(MobiDB)', 'database(Pfam)',
		'database(InterPro)', 'database(SUPFAM)',
		'database(PROSITE)', 'database(MUTAGENESIS)',
		#functional
		'feature(ACTIVE%20SITE)', 'feature(BINDING%20SITE)', 
		'feature(DNA%20BINDING)', 'feature(METAL%20BINDING)',
		'feature(NP%20BIND)', 'feature(SITE)',
		#topologies
		#INTRAMEM - Extent of a region located in a membrane without crossing it.
		'feature(TRANSMEMBRANE)', 'feature(INTRAMEMBRANE)', 
		'feature(TOPOLOGICAL%20DOMAIN)',
		#PTMs
		#CHAIN - Extent of a polypeptide chain in the mature protein.
		#PEPTIDE - Extent of a released active peptide.
		'feature(LIPIDATION)', 'feature(MODIFIED%20RESIDUE)',
		'feature(CROSS%20LINK)', 
		'feature(DISULFIDE%20BOND)', 'feature(GLYCOSYLATION)',
		'feature(INITIATOR%20METHIONINE)', 'feature(PEPTIDE)',
		'feature(SIGNAL)', 'feature(TRANSIT)',
		#I don't think we need these
		#'feature(PROPEPTIDE)', 'feature(CHAIN)', 
		#secondary structure
		'feature(BETA%20STRAND)', 'feature(HELIX)', 'feature(TURN)',
		#Domains: Aren't we getting those from pfam?
		#COMPOSITIONAL%20BIAS - Extent of a compositionally biased region. Example: /note="Glu-rich (acidic)"
		#MOTIF - Short (up to 20 amino acids) sequence motif of biological interest.
		#REGION - Extent of a region of interest in the sequence. Examples: /note="Zymogen activation region", /note="Possesses antibiotic activity"
		'feature(COILED%20COIL)', 'feature(COMPOSITIONAL%20BIAS)',
		'feature(MOTIF)', 'feature(REGION)', 'feature(ZINC%20FINGER)',
		#perhaps for the header of the file 
		'families'
    ]
    
    search_string = ('https://www.uniprot.org/uniprot/' +
        f'?query=reviewed:{reviewed}' +
        f'+AND+accession:{uniprot_id}' + 
        '+&format=tab&columns=id,entry%20name,protein%20name,' +
        ','.join(features)
        )
    
    req2 = urllib.request.Request(search_string)
    with urllib.request.urlopen(req2) as f:
       response2 = f.read()
    result = [i.split("\t") for i in response2.decode("utf-8").split("\n")]
    #retain the same terms for the column headers as were queried to reduce confusion
    result = pd.DataFrame(data=result[1:], columns=['Entry','Entry name']+features)
    
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
