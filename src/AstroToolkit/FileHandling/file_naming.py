def name_file(struct):
    if hasattr(struct, "subkind"):
        suffix = struct.subkind
    else:
        suffix = struct.kind

    if hasattr(struct, "survey"):
        if struct.survey:
            survey_str = f"{struct.survey}_"
        else:
            survey_str = ""
    else:
        survey_str = ""

    if struct.kind == "image":
        extension = ".fits"
    else:
        extension = ".csv"

    if hasattr(struct, "source"):
        if struct.source:
            fname = f"{struct.identifier}_{struct.source}_{survey_str}ATK{suffix}{extension}"
        elif struct.pos:
            fname = f"{struct.identifier}_{survey_str}ATK{suffix}{extension}"
    elif hasattr(struct, "sources"):
        if len(struct.identifiers) > 1:
            fname = f"{struct.identifiers[0]}_AndOtherSource(s)_ATK{suffix}{extension}"
        else:
            fname = f"{struct.identifiers[0]}_ATK{suffix}{extension}"

    return fname


def generate_plotname(struct, subkind=None):
    dataname = struct.dataname
    if subkind:
        dataname = dataname.replace(f"ATK{struct.kind}", subkind)

    if dataname.endswith(".csv"):
        struct.plotname = dataname[:-4] + ".html"
    elif dataname.endswith(".fits"):
        struct.plotname = dataname[:-5] + ".html"
