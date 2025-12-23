import ATK.Tools as ATK

SOURCE = 2552928187080872832
SOURCE = 587316166180416640

target = ATK.Target.from_id(SOURCE)

data = ATK.query("image", survey="panstarrs", target=target, size=150, band="r", overlays=["gaia", "galex"])

# data = ATK.query("image", survey="dss1", target=SOURCE, size=300, band="blue", overlays=["galex"])

if data.data:
    img = data.data[0]
    img.show()

    print("\n\n---------------------------------\n\n")

    data.show()
    data.plot(relative_axes=True)
    data.open()
