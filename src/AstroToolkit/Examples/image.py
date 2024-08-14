from AstroToolkit.Tools import query

# specify a Gaia source
source = 587316166180416640

""" 
Retrieve any available image and plot it.
No size or band is given, so these will be taken from the config. 
Since overlays is given as a list, only the magnitude listed for each survey in the config
will be overlayed as a detection (in this case, phot_g_mean_mag for gaia, and NUVmag for GALEX)
"""
figure = query(
    kind="image", source=source, survey="any", overlays=["gaia", "galex"]
).plot()

# show the image in the browser (and save to a static .html file)
figure.showplot()

# Now, give overlays as a dict containing all magnitudes to overlay.
figure_allmags = query(
    kind="image",
    source=source,
    survey="any",
    overlays={
        "gaia": ["phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag"],
        "galex": ["NUVmag", "FUVmag"],
    },
).plot()

# show the image in the browser (and save to a static .html file)
figure_allmags.showplot()

"""
Now, include a light curve survey as an overlay. For surveys with enough positional
precision, this can be used to trace the motion of the object through time.
We have also specified an image size, which will override the value in the config.
A different gaia source with a large proper motion has been used.
"""
figure_tracer = query(
    kind="image",
    source=2552928187080872832,
    survey="panstarrs",
    overlays=["crts"],
    size=60,
).plot()

figure_tracer.showplot()
