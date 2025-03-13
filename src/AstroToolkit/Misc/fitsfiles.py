from difflib import SequenceMatcher

from astropy.table import Table

newline = "\n"


def get_columns(filename, columns):
    if not filename.endswith(".fits"):
        filename += ".fits"
    if not isinstance(columns, list):
        columns = [columns]

    try:
        returned_data = Table.read(filename).to_pandas()
    except:
        raise Exception(f"Note: {filename} not found, or invalid format")

    data = []
    for column in columns:
        try:
            data.append(returned_data.loc[:, column].tolist())
        except:
            headers = returned_data.columns.values.tolist()
            similarity_arr = []
            for header in headers:
                similarity_arr.append(SequenceMatcher(None, column, header).ratio())

            similarity_arr, headers = zip(*sorted(zip(similarity_arr, headers)))
            raise Exception(
                f'No "{column}" column found.{newline}{newline}Possible similar columns found: {headers[-1]}, {headers[-2]}, {headers[-3]}'
            )

    return data
