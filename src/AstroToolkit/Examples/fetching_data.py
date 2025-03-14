from AstroToolkit.Tools import query

gaia_data = query(kind="data", source=587316166180416640, survey="gaia")

galex_data = query(kind="data", source=587316166180416640, survey="galex")
galex_data.showdata()

print("\n", galex_data.data["FUVmag"][0], galex_data.data["NUVmag"][0])
