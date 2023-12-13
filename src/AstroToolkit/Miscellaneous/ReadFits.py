from astropy.table import Table

def get_source_list(file_name):
    if not file_name.endswith('.fits'):
        file_name+='.fits'

    try:
        data=Table.read(file_name).to_pandas()
    except:
        raise Exception(f'File: {file_name} not found, or invalid format.')
    
    try:
        sources=data.loc[:,'source_id'].tolist()
    except:
        raise Exception('No source_id column found.')

    return sources

def get_pos_list(file_name):
    if not file_name.endswith('.fits'):
        file_name+='.fits'

    try:
        data=Table.read(file_name).to_pandas()
    except:
        raise Exception(f'File: {file_name} not found, or invalid format.')
    
    try:
        ra=data.loc[:,'ra'].tolist()
        dec=data.loc[:,'dec'].tolist()
    except:
        raise Exception('ra/dec columns not found.')
    
    pos_list=[list(x) for x in zip(ra,dec)]
    return pos_list