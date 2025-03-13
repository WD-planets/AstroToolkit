from AstroToolkit.Datapages import buttons, datapage, datatable
from AstroToolkit.Tools import query

# source = Hu Leo
source = 587316166180416640

# get image data and plot it
image = query(
    kind="image",
    survey="panstarrs",
    source=source,
    overlays=["gaia", "galex"],
    check_exists="datapage_image",
).plot()

# get hrd data and plot it
hrd = query(kind="hrd", sources=source, check_exists="datapage_hrd").plot()

# get spectrum data and plot it
spectrum = query(
    kind="spectrum", survey="sdss", source=source, check_exists="datapage_spectrum"
).plot()

# get sed data and plot it
sed = query(kind="sed", source=source, check_exists="datapage_sed").plot(
    spectrum_overlay=spectrum
)

# get lightcurve data [g,r,i] and plot it
lightcurves = query(
    kind="lightcurve", survey="ztf", source=source, check_exists="datapage_lightcurve"
).plot(colours=["green", "red", "blue"])

# plot a power spectrum (doesn't re-perform query as data already exists)
powspec = query(
    kind="lightcurve", survey="ztf", source=source, check_exists="datapage_lightcurve"
).plot(kind="powspec")

# get SIMBAD and Vizier buttons
buttons = buttons(source=source)

# get a metadata table with default parameters for various surveys
metadata = datatable(
    source=source,
    entries=[
        {
            "kind": "atk_defaults",
            "surveys": [
                "gaia",
                "galex",
                "panstarrs",
                "skymapper",
                "sdss",
                "wise",
                "twomass",
            ],
        }
    ],
)

# formats plots for use in grid
datapage = datapage(
    dimensions={"width": 6, "height": 6},
    panels=[
        {"name": "image", "figure": image, "width": 2, "height": 2},
        {"name": "hrd", "figure": hrd, "width": 2, "height": 2},
        {"name": "sed", "figure": sed, "width": 3, "height": 2},
        {"name": "buttons", "figure": buttons, "width": 2, "height": 2},
        {"name": "lightcurves", "figure": lightcurves, "width": 4, "height": 2},
        {"name": "powspec", "figure": powspec, "width": 3, "height": 2},
        {"name": "spectrum", "figure": spectrum, "width": 5, "height": 2},
        {"name": "metadata_table", "figure": metadata, "width": 7, "height": 2},
    ],
    layout=[
        ["image", "sed", "buttons"],
        ["hrd", "spectrum"],
        ["lightcurves", "powspec"],
        ["metadata_table"],
    ],
    grid_size=200,
)

datapage.showplot(f"{source}_datapage")
