import ATK.Tools as ATK

SOURCE = 587316166180416640
# SOURCE = 2552928187080872832

target = ATK.Target.from_id(SOURCE)

data = ATK.query("image", survey="panstarrs", target=SOURCE, size=120, band="g", overlays=["sdss"], disable_corrections=True)

# data = ATK.query("image", survey="dss1", target=SOURCE, size=120, band="blue", overlays=["galex"])

data.save("test_image.fits")
rec_data = ATK.read("test_image.fits")
rec_data.show()
rec_data.plot(cmap="false_colour")
rec_data.open()
