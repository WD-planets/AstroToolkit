import ATK.Tools as ATK

SOURCE = 2552928187080872832

data = ATK.query("image", survey="panstarrs", target=SOURCE, size=120, band="r")
# data.save("test_image.fits")
# data.show()

# data = ATK.read("test_image.fits")
# data.show()

data.plot(cmap="false_colour")
data.open()
