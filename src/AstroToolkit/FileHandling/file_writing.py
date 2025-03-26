import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

newline = "\n"

pd.set_option("future.no_silent_downcasting", True)


def writeHeader(struct, fname, metadata):
    with open(fname, "w+") as f:
        f.write("# atk_metadata:\n")
        for key, val in metadata.items():
            f.write(f"# {key}={val}{newline}")
        f.write("\n")


def CreateLocalData(struct, fname):
    # Allows you to save data if no source or pos is given
    try:
        ra, dec = struct.pos[0], struct.pos[1]
    except:
        ra, dec = None, None

    data_dict = {}
    for key in struct.data:
        data_dict[key] = struct.data[key]

    df = pd.DataFrame.from_dict(data_dict)
    df = df.fillna(value=np.nan).infer_objects(copy=False)
    table = Table.from_pandas(df)
    hdu = fits.table_to_hdu(table)

    hdr = fits.Header()
    header_info = {
        "atk_kind": struct.kind,
        "atk_subkind": struct.subkind,
        "atk_survey": struct.survey,
        "atk_catalogue": struct.catalogue,
        "atk_source": struct.source,
        "atk_pos_ra": ra,
        "atk_pos_dec": dec,
        "atk_identifier": struct.identifier,
        "atk_t": struct.trace,
    }
    for key, val in header_info.items():
        hdr[key] = val
    hdu_list = fits.HDUList([fits.PrimaryHDU(header=hdr), hdu])
    hdu_list.writeto(fname, overwrite=True)

    return True


def CreateLocalBulkdata(struct, fname):
    try:
        ra, dec = struct.pos[0], struct.pos[1]
    except:
        ra, dec = None, None

    table_hdus = []
    for survey in struct.data:
        data_dict = {}
        if struct.data[survey]:
            for key in struct.data[survey]:
                data_dict[key] = struct.data[survey][key]

        df = pd.DataFrame.from_dict(data_dict)
        df = df.fillna(value=np.nan).infer_objects(copy=False)
        table = Table.from_pandas(df)
        table_hdu = fits.table_to_hdu(table)
        table_hdu.header["atk_survey"] = survey
        table_hdus.append(table_hdu)

    hdr = fits.Header()
    header_info = {
        "atk_kind": struct.subkind,
        "atk_subkind": struct.subkind,
        "atk_survey": survey,
        "atk_source": struct.source,
        "atk_pos_ra": ra,
        "atk_pos_dec": dec,
        "atk_identifier": struct.identifier,
        "atk_t": struct.trace,
    }
    for key, val in header_info.items():
        hdr[key] = val
    hdu_list = fits.HDUList(fits.PrimaryHDU(header=hdr))
    for hdu in table_hdus:
        hdu_list.append(hdu)
    hdu_list.writeto(fname, overwrite=True)

    return True


def CreateLocalLightcurve(struct, fname):
    try:
        ra, dec = struct.pos[0], struct.pos[1]
    except:
        ra, dec = None, None

    table_hdus = []
    for band in struct.data:
        data_dict = {}
        band_label = band["band"]
        if band["mag"] is not None:
            for key in band:
                if key != "band":
                    data_dict[key] = band[key]
        df = pd.DataFrame.from_dict(data_dict)
        df.fillna(value=np.nan).infer_objects(copy=False)
        table = Table.from_pandas(df)
        table_hdu = fits.table_to_hdu(table)
        table_hdu.header["atk_band"] = band_label
        table_hdus.append(table_hdu)

    hdr = fits.Header()
    header_info = {
        "atk_kind": struct.kind,
        "atk_survey": struct.survey,
        "atk_source": struct.source,
        "atk_pos_ra": ra,
        "atk_pos_dec": dec,
        "atk_identifier": struct.identifier,
        "atk_t": struct.trace,
    }
    for key, val in header_info.items():
        hdr[key] = val
    hdu_list = fits.HDUList(fits.PrimaryHDU(header=hdr))
    for hdu in table_hdus:
        hdu_list.append(hdu)
    hdu_list.writeto(fname, overwrite=True)

    return True


def CreateLocalImage(struct, fname):
    try:
        ra, dec = struct.pos[0], struct.pos[1]
    except:
        ra, dec = None, None

    from astropy.io import fits

    image_hdu = fits.PrimaryHDU(struct.data["image_data"], header=struct.data["image_header"])
    hdu_list = fits.HDUList(image_hdu)

    header_info = {
        "atk_kind": struct.kind,
        "atk_survey": struct.survey,
        "atk_source": struct.source,
        "atk_pos_ra": ra,
        "atk_pos_dec": dec,
        "atk_identifier": struct.identifier,
        "atk_image_focus_ra": struct.data["image_focus"][0],
        "atk_image_focus_dec": struct.data["image_focus"][1],
        "atk_image_size": struct.data["size"],
        "atk_image_time_year": struct.data["image_time"][0],
        "atk_image_time_month": struct.data["image_time"][1],
        "atk_t": struct.trace,
    }
    for key, val in header_info.items():
        hdu_list[0].header[key] = val

    if "overlay" in struct.data and struct.data["overlay"]:
        overlay_data = struct.data["overlay"]
        surveys = list(dict.fromkeys([x["survey"] for x in overlay_data]))

        for survey in surveys:
            overlay_dict = {}
            survey_data = [x for x in overlay_data if x["survey"] == survey]
            keys = list(survey_data[0].keys())
            for key in keys:
                data = []
                for entry in survey_data:
                    data.append(entry[key])
                overlay_dict[key] = data

            df = pd.DataFrame.from_dict(overlay_dict)
            table = Table.from_pandas(df)
            hdu = fits.table_to_hdu(table)
            hdu_list.append(hdu)

    hdu_list.writeto(fname, overwrite=True)

    return True


def CreateLocalReddening(struct, fname):
    success = CreateLocalData(struct, fname)
    return success


def CreateLocalSed(struct, fname):
    try:
        ra, dec = struct.pos[0], struct.pos[1]
    except:
        ra, dec = None, None

    hdr = fits.Header()
    header_info = {
        "atk_kind": struct.kind,
        "atk_source": struct.source,
        "atk_pos_ra": ra,
        "atk_pos_dec": dec,
        "atk_identifier": struct.identifier,
        "atk_t": struct.trace,
    }
    for key, val in header_info.items():
        hdr[key] = val

    data_dict = {"survey": [], "mag_name": [], "wavelength": [], "flux": [], "flux_rel_err": []}
    for data_set in struct.data:
        for i, _ in enumerate(data_set["wavelength"]):
            data_dict["survey"].append(data_set["survey"])
            data_dict["mag_name"].append(data_set["mag_name"][i])
            data_dict["wavelength"].append(data_set["wavelength"][i])
            data_dict["flux"].append(data_set["flux"][i])
            data_dict["flux_rel_err"].append(data_set["flux_rel_err"][i])

    df = pd.DataFrame.from_dict(data_dict)
    df = df.fillna(value=np.nan).infer_objects(copy=False)
    table = Table.from_pandas(df)
    hdu = fits.table_to_hdu(table)

    hdu_list = fits.HDUList([fits.PrimaryHDU(header=hdr), hdu])
    hdu_list.writeto(fname, overwrite=True)

    return True


def CreateLocalSpectrum(struct, fname):
    try:
        ra, dec = struct.pos[0], struct.pos[1]
    except:
        ra, dec = None, None

    hdr = fits.Header()
    header_info = {
        "atk_kind": struct.kind,
        "atk_survey": struct.survey,
        "atk_source": struct.source,
        "atk_pos_ra": ra,
        "atk_pos_dec": dec,
        "atk_identifier": struct.identifier,
        "atk_t": struct.trace,
    }
    for key, val in header_info.items():
        hdr[key] = val

    df = pd.DataFrame.from_dict(struct.data)
    df = df.fillna(value=np.nan).infer_objects(copy=False)
    table = Table.from_pandas(df)
    hdu = fits.table_to_hdu(table)

    hdu_list = fits.HDUList([fits.PrimaryHDU(header=hdr), hdu])
    hdu_list.writeto(fname, overwrite=True)

    return True


def CreateLocalHrd(struct, fname):
    hdr = fits.Header()
    header_info = {"atk_kind": struct.kind, "atk_survey": struct.survey, "atk_t": struct.traces}
    for key, val in header_info.items():
        hdr[key] = val

    df = pd.DataFrame.from_dict(struct.data)
    df["sources"] = struct.sources
    df["identifiers"] = struct.identifiers
    df["position_ra"] = [x[0] for x in struct.positions]
    df["position_dec"] = [x[1] for x in struct.positions]
    df = df.fillna(value=np.nan).infer_objects(copy=False)
    table = Table.from_pandas(df)
    hdu = fits.table_to_hdu(table)

    hdu_list = fits.HDUList([fits.PrimaryHDU(header=hdr), hdu])
    hdu_list.writeto(fname, overwrite=True)

    return True


def generate_local_file(struct, name):
    if hasattr(struct, "subkind"):
        ftype = struct.subkind
    else:
        ftype = struct.kind

    if not struct.data:
        print("Note: No data to save, suggests that no data was found.")
        return None

    if not name:
        fname = struct.dataname
    else:
        fname = name

    if not struct.dataname:
        raise ValueError(".dataname attribute of structure is None.")

    success = globals()[f"CreateLocal{ftype.capitalize()}"](struct, fname)
    return success
