import astropy.units as u
import cairosvg
from astropy.coordinates import SkyCoord
from bokeh.io import export_png, export_svg, output_file
from utilities import format_plot

from ATK import query

vMa2 = SkyCoord(12.291, 5.389, unit="deg", frame="icrs")
data = query("vizier", targets=vMa2, survey="galex", path="vizier_1.fits.gz")
data.show()

vMa2 = 2552928187080872832
data = query("vizier", targets=vMa2, survey="galex", path="vizier_2.fits.gz")
data.show()

vMa2 = SkyCoord(12.291, 5.389, unit=u.deg, frame="icrs")
data = query(
    "image",
    targets=vMa2,
    survey="panstarrs",
    band="y",
    size=2 * u.arcmin,
    overlays=["galex"],
    disable_correction=True,
    path="image_1.fits.gz",
)
data.plot()
data.figure = format_plot(data.figure, 2, 2, True)
export_png(data.figure, filename="uncorrected_image.png")

vMa2 = 2552928187080872832
data = query("image", targets=vMa2, survey="panstarrs", band="y", size=2 * u.arcmin, overlays=["galex"], path="image_2.fits.gz")
data.show()
data.plot()
data.figure = format_plot(data.figure, 2, 2, True)
export_png(data.figure, filename="corrected_image.png")

ARSco = 6050296829033196032
data = query("lightcurve", targets=ARSco, survey="asassn", path="lightcurve.fits.gz")
data.show()
data.plot()
data.figure = format_plot(data.figure, 2, 1)
export_png(data.figure, filename="lightcurve.png")

pspec = data.apply("pspec", fmin=0, fmax=10, samples=100000, inplace=False)
pspec.show()
pspec.plot()
pspec.figure = format_plot(pspec.figure, 2, 2, True)
export_png(pspec.figure, filename="pspec.png")

fold = data.apply("fold", fmin=0, fmax=10, samples=100000, inplace=False)
fold.show()
fold.plot()
fold.figure = format_plot(fold.figure, 2, 1)
export_png(fold.figure, filename="folded_lightcurve.png")

HuLeo = 587316166180416640
spec = query("spectrum", targets=HuLeo, survey="sdss", path="spectrum.fits.gz")
spec.show()
spec.plot()
spec.figure = format_plot(spec.figure, 2, 1)
export_png(spec.figure, filename="spec.png")

fitted = spec.plot(fit=True, smooth=5.0)
fitted.figure = format_plot(fitted.figure, 2, 1)
export_png(fitted.figure, filename="fitted_spec.png")

sed = query("sed", targets=HuLeo, path="sed.fits.gz")
sed.plot()
sed.figure = format_plot(sed.figure, 2, 1)
export_png(sed.figure, filename="sed.png")

hrd = query("hrd", targets=HuLeo, path="hrd.fits.gz")
hrd.plot()
hrd.figure = format_plot(hrd.figure, 2, 2, True)
export_png(hrd.figure, filename="hrd.png")
