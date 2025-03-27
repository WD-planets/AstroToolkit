from ..Configuration.epochs import EpochStruct

epochs = EpochStruct().epoch_list


class Correct(object):
    """Corrects coordinates using pmra and pmdec in mas/yr"""

    def __init__(self, input_time, target_time, ra, dec, pmra, pmdec):
        self.input_time = input_time
        self.target_time = target_time
        self.ra = ra
        self.dec = dec
        self.pmra = pmra
        self.pmdec = pmdec

    # check if pmra or pmdec are nan
    @property
    def check_nans(self):
        import math

        if math.isnan(self.pmra) or math.isnan(self.pmdec):
            print("Note: could not retrieve object's pmra/pmdec, and so coordinates/radius were not corrected.")
            return False
        else:
            return True

    # return the gap in years/months between input_time and target_time
    @property
    def get_deltas(self):
        input_year, input_month = self.input_time
        target_year, target_month = self.target_time

        self.year_delta = target_year - input_year
        self.month_delta = target_month - input_month

    # perform correction
    @property
    def correction(self):
        import math

        self.ra += (
            (self.year_delta * self.pmra / 3600000 + self.month_delta * self.pmra / 43200000)
            * 1
            / math.cos(self.dec / 360 * 2 * math.pi)
        )
        self.dec += self.year_delta * self.pmdec / 3600000 + self.month_delta * self.pmdec / 43200000

        return [self.ra, self.dec]


class CorrectRadius(Correct):
    """Corrects a radius using pmra and pmdec in mas/yr"""

    def __init__(self, input_time, target_time, radius, pmra, pmdec):
        self.input_time = input_time
        self.target_time = target_time
        self.radius = radius
        self.pmra = pmra
        self.pmdec = pmdec

    # perform correction
    @property
    def correction(self):
        import math

        time_delta = abs(self.year_delta + self.month_delta / 12)
        self.radius += math.sqrt((self.pmra / 1000) ** 2 + (self.pmdec / 1000) ** 2) * time_delta

        return self.radius


"""Corrects coordinates for proper motion given in mas/yr"""


def correctpm(input_time, target_time, ra, dec, pmra, pmdec, check_success=False):
    input = Correct(input_time, target_time, ra, dec, pmra, pmdec)
    check = input.check_nans
    if check:
        input.get_deltas
        corrected_coords = input.correction
        if check_success:
            return corrected_coords, True
        else:
            return corrected_coords
    else:
        if check_success:
            return [ra, dec], False
        else:
            return [ra, dec]


"""Corrects a radius for proper motion given in mas/yr"""


def correctradius(source, input_time, target_time, radius, check_success=False):
    from ..Tools import query

    gaia_data = query(kind="data", source=source, survey="gaia", radius=3, level="internal").data
    if gaia_data:
        pmra, pmdec = gaia_data["pmra"][0], gaia_data["pmdec"][0]
    else:
        print("Note: radius correction failed. Some detections may be missing for high proper motion systems.")
        return radius

    input = CorrectRadius(input_time, target_time, radius, pmra, pmdec)
    check = input.check_nans
    if check:
        input.get_deltas
        corrected_radius = input.correction
        if check_success:
            return corrected_radius, True
        else:
            return corrected_radius
    else:
        if check_success:
            return radius, False
        else:
            return radius


"""Performs automatic proper motion correction between two surveys for a Gaia source"""


def autocorrect_pos(input_survey, target_survey, ra=None, dec=None, pmra=None, pmdec=None, check_success=False):
    input_time, target_time = epochs[input_survey], epochs[target_survey]
    return correctpm(input_time, target_time, ra, dec, pmra, pmdec, check_success=check_success)


def autocorrect_source(source, target_time=None, target_survey=None, check_success=False):
    from ..Tools import query

    if target_survey:
        target_time = epochs[target_survey]

    gaia_data = query(kind="data", survey="gaia", source=source, radius=3, level="internal").data
    ra, dec, pmra, pmdec = (gaia_data["ra"][0], gaia_data["dec"][0], gaia_data["pmra"][0], gaia_data["pmdec"][0])
    return correctpm(epochs["gaia"], target_time, ra, dec, pmra, pmdec, check_success=check_success)
