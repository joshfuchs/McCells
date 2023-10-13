import cellxgene_census as cell_census
import anndata as ad
import time
from pronto import Ontology
import pandas as pd


def get_data(organism, obs_val_filter, var_val_filter, col_names):
    '''
    This function queries data from an active cell census object.
    
    Assumes there is an active census object already open. Assumes you only want to query on cell metadata. 
    Gene metadata querying not currently supported.
    
    Parameters
    ----------
    organism : string
        string containing type of organism you want to query for, i.e. 'Homo sapiens'
        
    obs_val_filter : string
        string containing obs parameters values to filter on
        
    col_names : dictionary
        dictionary of obs or val values to return
        
    Returns
    -------
        AnnData object
    
    '''
    adata = cell_census.get_anndata(
        census=census,
        organism = organism,
        obs_value_filter = obs_val_filter,
        var_value_filter = var_val_filter,
        column_names= col_names,
        )
    return adata


#cl = Ontology.from_obo_library('cl.owl')

# load gene list
biomart = pd.read_csv('mart_export.txt')

coding_only = biomart[biomart['Gene type'] == 'protein_coding']

gene_list = coding_only['Gene stable ID'].to_list()

var_val_filter = '''feature_id in {}'''.format(gene_list)


# specify target cell IDs
#leaf_list_in_cl = ['CL:0002343', 'CL:3000001', 'CL:0000895', 'CL:0000900', 'CL:0002394', 'CL:0002399', 'CL:2000055', 'CL:0001050', 'CL:0001044', 'CL:0000807', 'CL:0000808', 'CL:0000794', 'CL:0000985', 'CL:0000987', 'CL:0000913', 'CL:0000905', 'CL:0000091', 'CL:0000903', 'CL:0000910', 'CL:0000938', 'CL:0000915', 'CL:0000899', 'CL:0000904', 'CL:0002396', 'CL:0000934', 'CL:0000940', 'CL:0000907', 'CL:0001057', 'CL:0011025', 'CL:0000939', 'CL:0001062', 'CL:0001058', 'CL:0001043', 'CL:0001049', 'CL:0001076', 'CL:0002057'] #updated 24 August 2023 for leukocyte as top level


#leaf_list_in_cl = ['CL:0002343', 'CL:3000001', 'CL:0000895', 'CL:0000900', 'CL:0000233', 'CL:0002394', 'CL:0002399', 'CL:2000055', 'CL:0001050', 'CL:0001044', 'CL:0000807', 'CL:0000894', 'CL:0000808', 'CL:0000794', 'CL:0000985', 'CL:0000987', 'CL:0000913', 'CL:0000905', 'CL:0000091', 'CL:0000903', 'CL:0000910', 'CL:0001029', 'CL:0002355', 'CL:0002045', 'CL:0000559', 'CL:0000938', 'CL:0000936', 'CL:0002425', 'CL:0000915', 'CL:0000899', 'CL:0000904', 'CL:0002396', 'CL:0000934', 'CL:0000940', 'CL:0000907', 'CL:0001057', 'CL:0011025', 'CL:0000939', 'CL:0001062', 'CL:0001058', 'CL:0001043', 'CL:0001049', 'CL:0001076', 'CL:0002057'] # from 14 Sep 2023 for hematopoietic cell (988) as top cell
leaf_list_in_cl = ['CL:0001044','CL:0000624','CL:0000895','CL:0000576','CL:0002057']

obs_val_filter = '''assay == "10x 3\' v3" and cell_type_ontology_term_id in {}'''.format(leaf_list_in_cl)

organism = "Homo sapiens"
col_names = {"obs": ["cell_type_ontology_term_id"]}

print('starting download with obs filter of')
print(obs_val_filter)

start = time.time()

with cell_census.open_soma(census_version="stable") as census:
    adata = get_data(organism, obs_val_filter, var_val_filter, col_names)
    adata.write(filename='13Oct_2int_3leaf',compression='gzip')


end = time.time()
print(end-start)

