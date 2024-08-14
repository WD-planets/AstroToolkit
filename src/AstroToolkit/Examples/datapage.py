from bokeh.layouts import column, layout, row
from bokeh.plotting import output_file, show

from AstroToolkit.Datapages import buttons, datatable, gridsetup
from AstroToolkit.Tools import query

# source = Hu Leo
source = 587316166180416640

# set grid size (scales size of datapage)
grid_size = 250

# get image data and plot it
image = query(
    kind="image", survey="any", source=source, overlays=["gaia", "galex"]
).plot()

# get hrd data and plot it
hrd = query(kind="hrd", sources=source).plot()

# get sed data and plot it
sed = query(kind="sed", source=source).plot(spectrum_overlay=True, survey="sdss")

# get spectrum data and plot it
spectrum = query(kind="spectrum", survey="sdss", source=source).plot()

# get lightcurve data [g,r,i] and plot it
lightcurves = query(kind="lightcurve", survey="ztf", source=source).plot(
    colours=["green", "red", "blue"]
)

# get SIMBAD and Vizier buttons
buttons = buttons(source=source, grid_size=grid_size)

# get a metadata table with default parameters for various surveys
metadata = datatable(
    source=source,
    selection={
        "gaia": "default",
        "galex": "default",
        "panstarrs": "default",
        "skymapper": "default",
        "sdss": "default",
        "wise": "default",
        "twomass": "default",
    },
)

# formats plots for use in grid
grid_plots = gridsetup(
    dimensions={"width": 6, "height": 6},
    plots=[
        {"name": "image", "figure": image, "width": 2, "height": 2},
        {"name": "hrd", "figure": hrd, "width": 2, "height": 2},
        {"name": "sed", "figure": sed, "width": 2, "height": 2},
        {"name": "lightcurves", "figure": lightcurves, "width": 2, "height": 1},
        {"name": "buttons", "figure": buttons, "width": 2, "height": 1},
        {"name": "spectrum", "figure": spectrum, "width": 4, "height": 2},
        {"name": "metadata_table", "figure": metadata, "width": 6, "height": 2},
    ],
    grid_size=grid_size,
)

# set up the final grid
datapage = layout(
    column(
        row(
            grid_plots["image"],
            grid_plots["hrd"],
            column(grid_plots["buttons"], grid_plots["lightcurves"]),
        ),
        row(grid_plots["sed"], grid_plots["spectrum"]),
        row(grid_plots["metadata_table"]),
    )
)

# give output file a name
output_file(f"{source}_datapage.html")

# show the datapage (also saves it)
show(datapage)
