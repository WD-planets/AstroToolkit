import ATK.Tools as ATK

SOURCE = 2552928187080872832

# data = ATK.query("image", survey="panstarrs", target=SOURCE, size=120, band="r")

data = ATK.query("image", survey="sdss", target=SOURCE, size=120, band="g")

# data.save("test_image.fits")
# data.show()

# data = ATK.read("test_image.fits")
# data.show()

data.show()
data.plot(relative_axes=False)
data.open()
