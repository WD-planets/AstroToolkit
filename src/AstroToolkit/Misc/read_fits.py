from astropy.table import Table
from difflib import SequenceMatcher

newline='\n'

def get_column(file_name,parameter):
    # add file extension if missing
    if not file_name.endswith('.fits'):
        file_name+='.fits'

    # read fits table
    try:
        data=Table.read(file_name).to_pandas()
    except:
        raise Exception(f'File: {file_name} not found, or invalid format.')
    
    # search for column with given parameter name
    try:
        values=data.loc[:,parameter].tolist()
        return values
    except:
        # print top 3 similar column names if no column found with given name
        columns=data.columns.values.tolist()
        similarity_arr=[]
        for col in columns:
            similarity_arr.append(SequenceMatcher(None,parameter,col).ratio())
        
        similarity_arr,columns=zip(*sorted(zip(similarity_arr,columns)))
        raise Exception(f'No "{parameter}" column found.{newline}{newline}Possible similar columns found: {columns[-1]}, {columns[-2]}, {columns[-3]}')

def get_source_list(file_name,parameter):
    # add file extension if missing
    if not file_name.endswith('.fits'):
        file_name+='.fits'

    # read fits table
    try:
        data=Table.read(file_name).to_pandas()
    except:
        raise Exception(f'File: {file_name} not found, or invalid format.')

    # search for source column using given name
    try:
        sources=data.loc[:,parameter].tolist()
        return sources
    except:
        # print top 3 similar column names if no column found in given name
        columns=data.columns.values.tolist()
        similarity_arr=[]
        for col in columns:
            similarity_arr.append(SequenceMatcher(None, parameter,col).ratio())

        similarity_arr,columns=zip(*sorted(zip(similarity_arr,columns)))
        raise Exception(f'No "{parameter}" column found.{newline}{newline}Possible similar columns found: {columns[-1]}, {columns[-2]}, {columns[-3]}')

def get_pos_list(file_name,parameters):
    if not file_name.endswith('.fits'):
        file_name+='.fits'

    try:
        data=Table.read(file_name).to_pandas()
    except:
        raise Exception(f'File: {file_name} not found, or invalid format.')

    try:
        ra=data.loc[:,parameters[0]].tolist()
        dec=data.loc[:,parameters[1]].tolist()
    except:
        columns=data.columns.values.tolist()
        ra_similarity_arr=[]
        dec_similarity_arr=[]
        for col in columns:
            ra_similarity_arr.append(SequenceMatcher(None, parameters[0],col).ratio())
            dec_similarity_arr.append(SequenceMatcher(None, parameters[1],col).ratio())

        ra_similarity_arr,ra_columns_sorted=zip(*sorted(zip(ra_similarity_arr,columns)))
        dec_similarity_arr,dec_columns_sorted=zip(*sorted(zip(dec_similarity_arr,columns)))
        raise Exception(f'No "{parameters}" columns found.{newline}{newline}Possible similar columns found:{newline}{ra_columns_sorted[-1]}, {ra_columns_sorted[-2]}, {ra_columns_sorted[-3]}{newline}{dec_columns_sorted[-1]},{dec_columns_sorted[-2]},{dec_columns_sorted[-3]}')
    
    pos_list=[list(x) for x in zip(ra,dec)]
    
    return pos_list