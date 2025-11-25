from ATK.Tools import query

SOURCE = 2552928187080872832

data = query(kind="vizier", target=SOURCE, survey="gaia")
# data = query(kind="vizier", target=SOURCE, survey="galex")

for key, val in data.__dict__.items():
    print(key, val)
