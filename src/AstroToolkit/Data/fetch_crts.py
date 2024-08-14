import os
import time
from datetime import datetime
from io import StringIO

import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from bs4 import BeautifulSoup as bs
from importlib_resources import files
from numpy import extract, genfromtxt

TIMER = files("AstroToolkit.Configuration").joinpath("CRTS_TIMER.txt")

if not os.path.isfile(str(TIMER)):
    with open(str(TIMER), "w") as timer:
        current_time = datetime.now().strftime("%H:%M:%S")
        timer.write(current_time)


# Adapted from code given by Keith, as CRTS interface is not well-suited to scripted queries
def get_CRTS_lightcurve(pos, radius):
    ra, dec = pos[0], pos[1]

    survey_radius = int(radius / 60)

    with open(str(TIMER), "r") as timer:
        prev_query_time = timer.readline().split(":")

    # calculate time delta in seconds between now and time of last query
    current_time = datetime.now().strftime("%H:%M:%S").split(":")

    # this should fix the case of e.g. 23.59.59 turning to 00.00.01 the next day, which would then report a huge change in time and would allow two queries to be performed in a short time.
    # since there is no differentiation between days, can have issue with it seeing same time 24 hrs apart, and therefore making query wait 15s, but it doesn't matter enough to warrant fixing
    if current_time[0] == 0 and prev_query_time[0] == 23:
        delta_seconds = 0
    else:
        current_time_seconds = (
            int(current_time[2])
            + 60 * int(current_time[1])
            + 3600 * int(current_time[0])
        )
        prev_query_time_seconds = (
            int(prev_query_time[2])
            + 60 * int(prev_query_time[1])
            + 3600 * int(prev_query_time[0])
        )

        delta_seconds = abs(current_time_seconds - prev_query_time_seconds)

    # limit frequency of queries
    if delta_seconds < 15:
        wait_time = 15 - delta_seconds
        print(
            f"CRTS queries are limited to one per 15s. Waiting {wait_time}s before performing query."
        )
        time.sleep(wait_time)

    requestcoords = SkyCoord(ra, dec, unit="deg", frame="icrs")
    # Get CRTS data   http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv_new.cgi
    # The old url  = http://nunuku.caltech.edu/cgi-bin/getcssconedb_release_img.cgi?RA=+str(ra)+"&Dec="+str(dec)+"&Rad="+str(radius)+"&DB=photcat&OUT=csv&SHORT=long&PLOT=no"

    # radius / 60 converts from arcseconds to arcmin
    url = f"http://nunuku.caltech.edu/cgi-bin/getcssconedb_priv_new.cgi?RADec={ra} {dec}&Rad={survey_radius}&OUT=csv&SHORT=short&DB=photocat"

    with open(str(TIMER), "w") as timer:
        current_time = datetime.now().strftime("%H:%M:%S")
        timer.write(current_time)

    try:
        cqresult = requests.get(url)
        cqp = bs(cqresult.content, "html.parser")
        filelink = [link.get("href") for link in cqp.find_all("a")]

        if filelink != []:
            filelink = filelink[0]
            y = requests.get(filelink)

            if y.status_code == 200:  # good transaction
                get_data = genfromtxt(
                    StringIO(y.content.decode(encoding="utf-8")),
                    delimiter=",",
                    names=True,
                    dtype=None,
                )

                # delete rows which are too far away - i.e. beyond search radius)
                return pd.DataFrame(
                    extract(
                        requestcoords.separation(
                            SkyCoord(
                                get_data["RA"],
                                get_data["Dec"],
                                unit="deg",
                                frame="icrs",
                            )
                        ).arcsecond
                        < radius,
                        get_data,
                    )
                )
            else:
                print("Note: Experiencing issues with CRTS")
                return None

        print("Note: CRTS lightcurvequery returned no data.")
        return None

    except:
        raise ConnectionError(
            "error in CRTS_light_curve ra=" + str(ra) + "  dec=" + str(dec)
        )
