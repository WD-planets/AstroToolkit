from AstroToolkit.Tools import query

image_data = query(
    kind="image",
    source=2552928187080872832,
    survey="panstarrs",
    overlays=["gaia", "galex"],
    size=120,
)

image_data.showdata()
image_data.plot()
image_data.showplot()
