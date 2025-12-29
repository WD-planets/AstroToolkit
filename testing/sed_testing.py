import ATK.Tools as ATK

# SOURCE = 2552928187080872832
SOURCE = 587316166180416640

target = ATK.Target.from_id(SOURCE)

data = ATK.query("sed", target=SOURCE)

data.show()
data.open()
