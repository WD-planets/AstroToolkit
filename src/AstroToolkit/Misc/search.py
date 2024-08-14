import webbrowser


class DbSearch(object):
    def __init__(self, radius, pos=None, source=None):
        self.pos = pos
        self.source = source
        self.radius = radius

    def get_params(self):
        from ..Data.dataquery import SurveyInfo
        from ..Tools import query
        from .pmcorrection import correctradius

        survey_times = SurveyInfo().times

        if self.source:
            self.radius = correctradius(
                source=self.source,
                input_time=survey_times["gaia"],
                target_time=[1990, 0],
                radius=self.radius,
            )

            gaia_data = query(
                kind="data", survey="gaia", source=self.source, level="internal"
            ).data
            if gaia_data:
                self.pos = [gaia_data["ra"][0], gaia_data["dec"][0]]
            else:
                raise Exception("Gaia source not found.")


class SimbadQuery(DbSearch):
    def do_search(self):
        url = f"https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={self.pos[0]}+{self.pos[1]}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius={self.radius}&Radius.unit=arcsec&submit=submit+query&CoordList="
        webbrowser.open_new_tab(url)


class VizierQuery(DbSearch):
    def do_search(self):
        url = f"https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={self.pos[0]}+{self.pos[1]}&-c.rs={self.radius}&-out.add=_r&-sort=_r&-out.max=$4"
        webbrowser.open_new_tab(url)


class DatapageButtons(DbSearch):
    def get_urls(self):
        return (
            f"https://simbad.cds.unistra.fr/simbad/sim-coo?Coord={self.pos[0]}+{self.pos[1]}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius={self.radius}&Radius.unit=arcsec&submit=submit+query&CoordList=",
            f"https://vizier.cds.unistra.fr/viz-bin/VizieR-4?-c={self.pos[0]}+{self.pos[1]}&-c.rs={self.radius}&-out.add=_r&-sort=_r&-out.max=$4",
        )


def do_search(kind, radius, pos=None, source=None):
    query_object = globals()[f"{kind.capitalize()}Query"](
        radius=radius, pos=pos, source=source
    )
    query_object.get_params()
    query_object.do_search()
