import ATK.Tools as ATK

# SOURCE = 2552928187080872832
SOURCE = 587316166180416640

target = ATK.Target.from_id(SOURCE)

data = ATK.query("sed", target=SOURCE)
data.save("test_sed.fits")
rec_data = ATK.read("test_sed.fits")
rec_data.show()
rec_data.open()
