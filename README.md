# AstroToolkit

AstroToolkit (ATK) is a set of useful tools for fetching, plotting, and analysing astronomical data.

## Table of Contents

1. [Installation](#installation)
2. [Introduction](#introduction)
3. [Configuration](#configuration)
4. [Tools](#tools)

### ATK Tools

- 1 [Raw Data Tools](#1-raw-data-tools)
	- 1.1 [dataquery](#11-dataquery)
	- 1.2 [photquery](#12-photquery)
	- 1.3 [bulkphotquery](#13-bulkphotquery)
	- 1.4 [getreddening](#14-getreddening)
- 2 [Imaging Tools](#2-imaging-tools)
	- 2.1 [imagequery](#21-imagequery)
	- 2.2 [plotimage](#22-plotimage)
- 3 [HRD Tools](#3-hrd-tools)
	- 3.1 [plothrd](#31-plothrd)
- 4 [Lightcurve Tools](#4-lightcurve-tools)
	- 4.1 [lightcurvequery](#41-lightcurvequery)
	- 4.2 [plotlightcurve](#42-plotlightcurve)
- 5 [SED Tools](#5-sed-tools)
	- 5.1 [sedquery](#51-sedquery)
	- 5.2 [plotsed](#52-plotsed)
- 6 [Spectrum Tools](#6-spectrum-tools)
	- 6.1 [spectrumquery](#61-spectrumquery)
	- 6.2 [plotspectrum](#62-plotspectrum)
- 7 [Timeseries Tools](#7-timeseries-tools)
	- 7.1 [plotpowspec](#71-plotpowspec)
	- 7.2 [tsanalysis](#72-tsanalysis)
- 8 [Datapage Tools](#8-datapage-tools)
	- 8.1 [gridsetup](#81-gridsetup)
	- 8.2 [getbuttons](#82-getbuttons)
	- 8.3 [getmdtable](#83-getmdtable)
- 9 [File Handling Tools](#9-file-handling-tools)
	- 9.1 [savedata](#91-savedata)
	- 9.2 [readdata](#92-readdata)
	- 9.3 [export](#93-export)
- 10 [Miscellaneous Tools](#10-miscellaneous-tools)
	- 10.1 [showplot](#101-showplot)
	- 10.2 [saveplot](#102-saveplot)
	- 10.3 [correctpm](#103-correctpm)
	- 10.4 [getdistance](#104-getdistance)
	- 10.5 [convfromdeg](#105-convfromdeg)
	- 10.6 [convtodeg](#106-convtodeg)
	- 10.7 [getcolumn](#107-getcolumn)
	- 10.8 [getsources](#108-getsources)
	- 10.9 [getpositions](#109-getpositions)

<br>

## Installation<a id="installation"></a>

The package can be installed as with any other package, e.g. using pip:
```
pip install AstroToolkit
```

Once the package has been installed, you should navigate to the package install location. This should be located in: **.../Lib/site-packages/AstroToolkit**, where **...** is your python install location. If you wish to find this, you can use the following command from the terminal:

```
python -c "from AstroToolkit.Tools import getpath; getpath()"
```

In this directory, run either buildwin.bat (Windows) or build.sh (Linux) depending on your operating system.

**NOTE:** See _README.txt_ (Linux) or _README - Windows.txt_ (Windows) in the above directory for any additional dependencies.

<br>

## Introduction<a id="introduction"></a>

ATK uses Bokeh as its primary plotting library. The official documentation can be found at [https://bokeh.org/](https://bokeh.org/). A key property of Bokeh plots is that they can be saved as static .html files, which can then be shared/accessed while retaining all interactivity.

Fundamentally, there are two types of tools in ATK:
1. Fetching, used to obtain data from a given survey
2. Plotting, used to plot the data acquired from the above

In all fetching tools, there are two possible ways to target your system of interest:
1. pos = [ra,dec] in degrees
2. source = Gaia Source ID

Where possible, it is usually best to use a source as input, as this enables a key feature of ATK: **Proper Motion Correction**.

A good example of this is in imaging queries. If a 'pos' is used as input, the result will simply be the image data returned by the chosen imaging survey at those exact coordinates. However, this may not be ideal in the case of an object with a large proper motion. If a source is used instead, the data returned will have accounted for this, resulting in the image being centered on the system in question. This concept is used throughout ATK when matching data from different surveys to a given system.

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/source_vs_pos.png?raw=true)

<br>

All ATK tools can be imported using:
```
from AstroToolkit.Tools import [tool name]
```

<br>

## Configuration<a id="configuration"></a>
Most default arguments within ATK can be modified through the use of a config file. This can be edited through use of the **editconfig** tool.

**Usage:**

```
editconfig(options)
```

where options is a dictionary containing key : value pairs to set in the config. The list of accepted keys and their default values are shown below:
- enable_notifications = False
	- If True, notifications denoting basic information as each ATK tool is executed will be shown in the terminal. Can be useful for tracking the flow of data.
- dataquery_radius = 3
	- This sets the default radius (in arcseconds) to use in the [dataquery](#11-dataquery) tool.
- photquery_radius = 3
	- This sets the default radius (in arcseconds) to use in the [photquery](#12-photquery) tool.
- bulkphotquery_radius = 3
	- This sets the default radius (in arcseconds) to use in the [bulkphotquery](#13-bulkphotquery) tool.
- imagequery_size = 30
	- This sets the default radius (in arcseconds) to use in the [imagequery](#21-imagequery) tool.
- imagequery_overlays = gaia
	- This sets the default radius (in arcseconds) to use in the [imagequery](#21-imagequery) tool.
- lightcurvequery_radius = 3
	- This sets the default radius (in arcseconds) to use in the [lightcurvequery](#31-lightcurvequery) tool.
- atlas_username = None
	- This sets the default username to use for ATLAS lightcurve queries in the [lightcurvequery](#31-lightcurvequery) tool.
- atlas_password = None
	- This sets the default password to use for ATLAS lightcurve queries in the [lightcurvequery](#31-lightcurvequery) tool.
- sed_radius = 3
	- This sets the default radius to use in the [sedquery](#41-sedquery) tool.
- spectrum_radius = 3
	- This sets the default radius to use in the [spectrumquery](#51-spectrumquery) tool.
- grid_size = 250
	- This sets the default grid size to use in the [gridsetup](#71-gridsetup) tool.
- button_simbad_radius = 3
	- This sets the default radius to use for the SIMBAD button in the [getbuttons](#72-getbuttons) tool.
- button_vizier_radius = 3
	- This sets the default radius to use for the Vizier button in the [getbuttons](#72-getbuttons) tool.
- plot_size = 400
	- This sets the default plot size for all plotting tools.
- readfits_sourcename = source_id
	- This sets the default column name to use in the [getsources](#108-getsources)

**Example:**

To set the default username and password to use for ATLAS queries:

```
from AstroToolkit.Tools import editconfig

editconfig({'atlas_username':'USERNAME','atlas_password':'PASSWORD'})
```

<br>

To see the current value of a config parameter, you can use the **getconfigvalue** tool. This will return the parameter's value, and print it to the terminal.

**Usage:**

```
getconfigvalue(parameter)
```

where:

```
parameter = string, name of config parameter from above list
```

**Returns**: value of this parameter in the config.

<br><br>

## Tools

In this section, the available tools will be outlined. Note that if a parameter is listed as having a default parameter CONFIG, this means that this parameter is taken from the config as listed above. These parameters can still be passed to the tool, in which case the config value will be ignored.

**Note:** most data returned from fetching/plotting tools takes the form of a dictionary. This dictionary contains the returned data, as well as basic information such as the pos/source used to acquire the data. This format is then used by other ATK functions (such as file showing/saving/reading).

**Note:** All plots created by ATK can have their legends toggled by double clicking the plot, and individual data can be hidden by clicking them in the legend.

<br>

**In all data returned by querying tools, the value of the 'data' key will be set to None if no data was returned. Likewise, in all plotting tools, the value of the 'plot' key will be set to None.**

<br>

### 1. Raw Data Tools<a id="1-raw-data-tools"></a>
### 1.1. dataquery<a id="11-dataquery"></a>

Returns all available data for a given survey (e.g. magnitudes, positions, etc.).

<br>

**Supported surveys:** gaia, galex, rosat, panstarrs, skymapper, sdss, twomass, wise, erosita

<br>

**Usage:**

```
dataquery(survey,pos=None,source=None,radius=CONFIG)
```

where:

```
survey = str, name of a supported survey
```
```
pos = list, [ra,dec]
```
```
source = int/str, Gaia source_id
```
```
radius = int/float, radius of query
```

**Returns:** dict

```
{
'survey' : str, survey of data
'type' : 'data'
'source' : int/str, source used to get data (None if a pos was used)
'pos' : [ra,dec], pos used to get data (None if a source was used)
'data' : dict, the returned data
}
```

<br>

**Example:**

To retrieve the parallax of a system through Gaia, and its WISE data:

```
from AstroToolkit.Tools import dataquery

source = 6050296829033196032

parallax = dataquery(survey='gaia',source=source)['data']['parallax']
wise_data = dataquery(survey='wise',source=source)['data']
```

<br><br>

### 1.2. photquery<a id="12-photquery"></a>

Returns data from a given survey, with columns filtered to only include photometry and other basic information.

<br>

**Supported surveys:** gaia,panstarrs,skymapper,galex,rosat,sdss,wise,twomass

<br>

**Usage:**

```
photquery(survey,pos=None,source=None,radius=CONFIG)
```

where:

```
survey = str, name of a supported survey
```
```
pos = list, [ra,dec]
```
```
source = int/str, Gaia source_id
```
```
radius = int/float, radius of query
```

**Returns:** dict

```
{
'survey' : str, survey of data
'type' : 'data'
'source' : int/str, source used to get data (None if a pos was used)
'pos' : [ra,dec], pos used to get data (None if a source was used)
'data' : dict, the returned data
}
```

<br>

**Example:**

To retrieve the 2MASS photometry for an object:
```
from AstroToolkit.Tools import photquery

data=photquery(survey='twomass',source=6050296829033196032)['data']
```

<br>

### 1.3. bulkphotquery<a id="13-bulkphotquery"></a>

Returns available photometry from all surveys supported by [photquery](#12-photquery).

<br>

**Usage:**

```
bulkphotquery(pos=None,source=None,radius=CONFIG)
```

where:

```
pos = list, [ra,dec]
```
```
source = int/str, Gaia source_id
```
```
radius = int/float, radius of query
```

**Returns:** dict

```
{
'type' : 'bulkphot'
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'data' : {
         'gaia' : dict, returned data (or None if no data returned)
         'galex' : dict, returned data (or None if no data returned)

          etc. for each survey in surveys supported by photquery
         }
}
```

**Note:** The value of each survey will be set to None if no data was returned. E.g. if no galex data was returned, the value of the ['galex'] key would be None.

<br>

**Example:**

To retrieve the gaia and galex data for an object:

```
from AstroToolkit.Tools import bulkphotquery

bulk_phot=bulkphotquery(source=6050296829033196032)['data']

gaia_data=bulk_phot['gaia']
galex_data=bulk_phot['galex']
```

<br>

### 1.4. getreddening<a id="14-getreddening"></a>

Returns the reddening of a Gaia source.

Currently supported reddening surveys:
- STILISM (https://stilism.obspm.fr/)

<br>

**Usage:**

```
getreddening(source)
```

where:
```
source = int/str, Gaia source_id
```

<br>

**Returns:** dict

```
{
'type' : 'reddening'
'source' : int/str, source used to get data (None if a pos was used)
'data' : {
         'dist' : distance to source (1/parallax)
         'red_dist' : actual distance to which survey's reddening estimate refers
         'red_dist_err' : error on red_dist
         'red' : survey's reddening estimate
         'red_upper' : reddening upper limit
         'red_lower' : reddening lower limit
         }
}
```

<br><br>

### 2. Imaging Tools<a id="2-imaging-tools"></a>

These functions retrieve and plotimages from supported surveys.

**Note:** When using a 'pos' as input, some detections can be missing for high proper motion objects. When instead using a source as input, this is no longer a problem as the detection search radius is increased to account for this proper motion.

<br>

### 2.1. imagequery<a id="21-imagequery"></a>

Retrieves an image from a given survey. Overlays  can be used to overlay detections from different surveys. There are two types of overlay: detections and tracers. Detections are proper motion-corrected circles that scale with the magnitude of the detection, and show that there is data available for a given object in a given survey. 

Tracers (which are detections created from lightcurve surveys) can be used both to show that lightcurve data is available, but also to trace an object through time.

**Note:** tracing an object through time works best with ZTF and CRTS overlays, as these give data with the necessary coordinate precision. Other tracer overlays are unlikely to work for this, and will hence only be useful to see where lightcurve data is available.

**Note:** In the case of ATLAS overlays, the radius used is very small due to the time taken to perform ATLAS queries. Only data close to the focus of the image (i.e. at the location of the target object) will be shown. 

**Example:**

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/tracer_example.png?raw=true)

<br>

Supported surveys:
- panstarrs, supported bands = g, r, i, z, y
- skymapper, supported bands = g, r, i, z, u, v
- dss, supported bands = g

**Note:** Can also use 'any' to perform an image query according to the hierarchy: panstarrs > skymapper > dss

<br>

Supported overlays: gaia, galex_nuv, galex_fuv, rosat,sdss, twomass, wise, ztf, erosita ,atlas, gaia_lc, asassn, crts

**Note:** Can also use 'all', which will enable overlays for all supported surveys.

<br>

**Usage:**

```
imagequery(survey,pos=None,source=None,size=CONFIG,band='g',overlays=CONFIG)
```

where:
```
survey = str, name of a supported survey
```
```
pos = list, [ra,dec]
```
```
source = int/str, Gaia source_id
```
```
size = int/float, size of image in arcsec
```
```
band = str, string containing the required bands (e.g. for all panstarrs bands, use band='grizy')
```
```
overlays = str, required detection overlays (e.g. for gaia and wise detections, use overlays='gaia,wise') 
```

**Returns:** dict

```
{
'type' : 'image'
'survey' : str, image survey
'source' : int/str, source used to get image (None if a pos was used)
'pos' : [ra,dec], pos used to get image (None if a source was used)
'data' : {
         'image_data' : array, image data
         'image_header' : astropy header, image header
         'location' : [ra,dec], actual location of the image
         'size' : int/float, image size in arcsec
         'image_time' : [year,month], image time
         'wcs' : astropy wcs object of image
         'overlay' : list of overlay entries
         }
}
```

**NOTE:**  overlays are stored as a list of individual detections in the format:

```
{
'survey' : str, survey of detection
'position: [ra,dec], coordinates of detection
'radius' : float, radius of detection
'corrected' : bool, whether or not the detection has been corrected for proper motion
'mag' : str, name of the magnitude (column heaader) from a given survey
'marker' : 'circle' or 'cross', detection symbol to overlay. Circles are scaled with radius, crosses are not (e.g. for surveys without a magnitude to scale by)
}
```

<br>

### 2.2. plotimage<a id="22-plotimage"></a>

Plots images in format returned by [imagequery](#21-imagequery).

<br>

**Usage:**

```
plotimage(data)
```

where:
```
data = dict in format returned by imagequery
```

**Returns:** dict

```
{
'type' : 'image'
'survey' : str, survey of image
'source' : int/str, source used to get image (None if a pos was used)
'pos' : [ra,dec], pos used to get image (None if a source was used)
'plot' : bokeh figure object, the actual plot
'ATKfilename' : the default filename given by ATK
}
```

**Note:** if a plot is not returned (e.g. if no data was supplied to the plotting tool), the 'plot' key will be set to None.

<br>

**Example:**

To retrieve and plot an image:
```
from AstroToolkit.Tools import imagequery,plotimage,showplot

image=imagequery(survey='any',source=6050296829033196032,overlays='gaia')
plot=plotimage(image)
showplot(plot)
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/imaging_example.png?raw=true)

<br><br>

### 3. HRD Tools<a id="3-hrd-tools"></a>
### 3.1. plothrd<a id="31-plothrd"></a>

Returns a HRD with a source or list of sources overlayed over a Gaia 100pc background sample.

<br>

**Usage:**

```
plothrd(source=None,sources=None)
```

where:

```
source = int/str, Gaia source_id
```
```
sources = list of sources
```

**Returns:** dict

```
{
'type' : 'hrd'
'source' : int/str, source overlayed in HRD
'sources' : list, sources overlayed in HRD if multiple were given
'plot' : bokeh figure object, the actual plot
'ATKfilename' : the default filename given by ATK
}
```

<br>

**Example:**
To retrieve a HRD with a single source overlayed:
```
from AstroToolkit.Tools import plothrd,showplot

plot=plothrd(source=6050296829033196032)
showplot(plot)
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/hrd_example.png?raw=true)

<br><br>

### 4. Lightcurve Tools<a id="4-lightcurve-tools"></a>
### 4.1. lightcurvequery<a id="41-lightcurvequery"></a>

Returns lightcurve data for a given survey.

<br>

Supported surveys:
- ZTF (ztf) - g, r, i
- ATLAS (atlas) - o, c, i
- ASAS-SN (asassn) - g, v
- Gaia (gaia) - g, bp, rp
- CRTS (crts) - v

**Note:** CRTS' data server is not designed for scripted queries. As a result, queries are limited in the toolkit to one query per 15 seconds. This cooldown is managed using the CRTS_TIMER.txt file in the Settings directory of the package files, and hence this should not be edited.

<br>

```
lightcurvequery(survey,pos=None,source=None,radius=CONFIG,username=CONFIG,password=CONFIG,sigmaclip=None)
```

where:
```
survey = str, name of a supported survey
```
```
pos = list, [ra,dec]
```
```
source = int/str, Gaia source_id
```
```
radius = int/float, radius of lightcurve query
```
```
username = str, ATLAS username, hence only used in ATLAS queries
```
```
password = str, ATLAS password, hence only used in ATLAS queries
```
```
sigmaclip = int, performs sigma clipping on the data to this number of standard deviations
```

**Returns:** list of dicts

list of lightcurve data dictionaries with an entry for each band in that survey (see below). Each entry has the format:

```
{
'type' : 'lightcurve'
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'survey' : str, survey of data
'band' : str, band of lightcurve data
'data' : {
         'ra' : list of returned ra values
         'dec' : list of returned dec values
         'hjd'/'mjd' : list of returned hjd/mjd values, with the minimum returned
          value subtracted from all values (i.e. is just a measure of days from the 
          first observation)
         'hjd_ori'/'mjd_ori' : list of returned hjd/mjd values, unedited
         'mag' : list of returned magnitude values
         'mag_err' : list of returned magnitude error values
         }
}
```

<br>

### 4.2. plotlightcurve<a id="42-plotlightcurve"></a>

Plots lightcurves in the format returned by [lightcurvequery](#41-lightcurvequery).

<br>

**Usage:**

```
plotlightcurve(data,colour='black')
```

where:
```
data = dict if passing a single lightcurve, list of dicts if passing multiple lightcurves
```
```
colour = str, name of a supported colour. Only used when passing a single lightcurve
```
```
colours = list of strings denoting supported colours, e.g. ['green','red','blue']. Only used when passing multiple lightcurves
```

**Returns:** dict

```
{
'type' : 'lightcurve'
'survey' : str, survey of lightcurve
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'plot' : bokeh figure object, the actual plot
'ATKfilename' : the default filename given by ATK
}
```

<br>

**Example**:
To retrieve and plot lightcurves from ZTF:
```
from AstroToolkit.Tools import lightcurvequery,plotlightcurve,showplot

lightcurves=lightcurvequery(survey='ztf',source=6050296829033196032)
showplot(plotlightcurve(lightcurves,colours=['green','red','blue']))
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/lightcurve_example.png?raw=true)

**Note:** Each band will then be toggleable using the legend.

<br><br>

### 5. SED Tools<a id="5-sed-tools"></a>
### 5.1. sedquery<a id="51-sedquery"></a>

Queries all supported photometry surveys and returns SED data.

<br>

**Usage:**

```
sedquery(pos=None,source=None,radius=CONFIG)
```

where:
```
pos = list, [ra,dec]
```
```
source = int/str, Gaia source_id
```
```
radius = int/float, radius of data query
```

**Returns:** dict

```
{
'type' : 'sed'
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'data': list of entries with each entry taking the form of a dict:
	{
    'survey' : str, survey of data point
    'wavelength' : filter wavelength of data point
    'flux' : flux through filter
    'rel_err' : relative error on flux
    }
}
```

<br>

### 5.2. plotsed<a id="52-plotsed"></a>

Plots SEDs in the format returned by [sedquery](#51-sedquery).

<br>

**Usage:**

```
plotsed(data)
```

where:

```
data = dict in format returned by sedquery
```

**Returns:** dict

```
{
'type' : 'sed'
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'plot' : bokeh figure object, the actual plot
'ATKfilename' : the default filename given by ATK
}
```

<br>

**Example:** To retrieve and plot an SED:
```
from AstroToolkit.Tools import sedquery,plotsed,showplot

data=sedquery(source=6050296829033196032)
showplot(plotsed(data))
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/sed_example.png?raw=true)

<br><br>

### 6. Spectrum Tools<a id="6-spectrum-tools"></a>
### 6.1. spectrumquery<a id="61-spectrumquery"></a>

Returns spectrum data from a given survey.

<br>

**Usage:**

```
spectrumquery(survey=None,pos=None,source=None,radius=CONFIG)
```

where:
```
survey = str, name of a supported survey
```
```
pos = list, [ra,dec]
```
```
source = int/str, Gaia source_id
```
```
radius = int/float, radius of data query
```

**Returns:** dict

```
{
'type' : 'spectra'
'survey' : str, survey of detection
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'data': {
        'wavelength' : list of wavelength values
        'flux' : list of flux values
        }
}
```

<br>

### 6.2. plotspectrum<a id="62-plotspectrum"></a>

Plots spectra in the format returned by [spectrumquery](#61-spectrumquery).

<br>

**Usage:**

```
plotspectrum(data)
```

where:

```
data = dict in format returned by spectrumquery
```

**Returns:** dict

```
{
'type' : 'spectrum'
'survey' : str, survey of spectrum
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'plot' : bokeh figure object, the actual plot
'ATKfilename' : the default filename given by ATK
}
```

<br>

**Example:**
To retrieve and plot an SDSS spectrum:
```
from AstroToolkit.Tools import spectrumquery,plotspectrum,showplot

data=spectrumquery(survey='sdss',source=587316166180416640)
plot=plotspectrum(data)
showplot(plot)
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/spectrum_example.png?raw=true)

<br><br>

### 7. Timeseries Tools<a id="7-timeseries-tools"></a>
### 7.1 plotpowspec<a id="71-plotpowspec"></a>

Plots a power spectrum from lightcurve data in the format returned by [lightcurvequery](#41-lightcurvequery)

<br>

**Usage:**

```
plotpowspec(data)
```

where:
```
data = dict or list in format returned by lightcurvequery
```

**Returns:** dict

```
{
'type' : 'powspec'
'survey' : str, survey of lightcurve data given to the powspec tool
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'plot' : bokeh figure object, the actual plot
'ATKfilename' : the default filename given by ATK
}
```

<br>

**Example:**
To retrieve lightcurve data from ZTF and plot a power spectrum:
```
from AstroToolkit.Tools import lightcurvequery,plotpowspec,showplot

data=lightcurvequery(survey='ztf',source=6050296829033196032)
plot=plotpowspec(data)
showplot(plot)
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/powspec_example.png?raw=true)

<br><br>

### 7.2 tsanalysis<a id="72-tsanalysis"></a>

Runs period analysis on lightcurve data in format returned by [lightcurvequery](#lightcurvequery)

<br>

**Usage:**

```
tsanalysis(data)
```

where:

```
data = dict or list in format returned by lightcurvequery
```

<br>

**Example:**
To retrieve lightcurve data from ZTF and perform timeseries analysis:
```
from AstroToolkit.Tools import lightcurvequery,tsanalysis

data=lightcurvequery(survey='ztf',source=6050296829033196032)
plot=tsanalysis(data)
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/tsanalysis_example.png?raw=true)

<br><br>

## 8. Datapage Tools<a id="8-datapage-tools"></a>

These functions are used to create custom datapages from any plots/data supported by AstroToolkit.

**NOTE:** An example of datapage creation can be found within the packages 'Examples' folder, named 'datapage_creation.py' (within the …/Lib/site-packages/AstroToolkit from earlier). This can be imported from a python terminal using from AstroToolkit.Examples import datapage_creation.

Below is the datapage produced by this example. Once the script to create these has been written, they can be a very powerful way to quickly retrieve and show a wide range of information on a system, with the below example easily being generated in under 30s.

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/datapage_example.png?raw=true)

<br>

### 8.1. gridsetup<a id="81-gridsetup"></a>

Formats plots for datapage creation by sizing them to a grid of given dimensions.

<br>

**Usage:**

```
getgrid(dimensions,plots,grid_size=CONFIG)
```

where:
```
dimensions = dict, dict with keys 'width' and 'height' that define the dimensions of the datapage in grid units (e.g. for a datapage that is 6 units wide and 3 units tall, dimensions={'width':6,'height':3})
```
```
plots = list of plots, with each entry having the format:
        {
        'name' : the label to give to the plot
        'figure' : the plot dictionary as returned by plotting tools
        'width' : the width in grid units of this element
        'height' : the height in grid units of this element
        }
```
```
grid_size  = int, size of each square of the grid to which all plots are scaled.
```

**Returns:** dict

```
{
'name 1' : plot 1
'name 2' : plot 2

etc. for each plot given to the tool
}
```

where the keys "name ..." are the names given to the plots above, and the values 'plot ...' are the plots to which these names refer.

**NOTE:** Again, see the datapage_creation example as noted [above](#8-datapage-tools) for an example.

<br>

### 8.2. getbuttons<a id="82-getbuttons"></a>

Returns a Bokeh figure containing SIMBAD and Vizier buttons for use in datapages.

<br>

**Usage:**

```
getinfobuttons(grid_size,source=None,pos=None,simbad_radius=CONFIG,vizier_radius=CONFIG)
```

where:
```
grid_size = int, size of the grid to which the buttons are scaled.
```
```
pos = list, [ra,dec] in degrees
```
```
source = int/str, Gaia source_id
```
```
simbad_radius = int, radius to use in SIMBAD queries
```
```
vizier_radius = int, radius to use in Vizier queries
```

**Returns:** dict

```
{
'type' : 'buttons'
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'plot' : bokeh figure object, the actual plot element
}
```

**NOTE:** Again, see the datapage_creation example as noted [above](#8-datapage-tools) for an example.

<br>

### 8.3. getmdtable<a id="83-getmdtable"></a>

Creates a table of metadata table using data from supported surveys and/or custom data.

<br>

**Usage:**

```
getmdtable(metadata,pos=None,source=None)
```

where:

```
metadata = dict, dictionary of metadata in accepted format (see below)
```
```
pos = list, [ra,dec] in degrees
```
```
source = int/str, Gaia source_id
```

<br>

The expected metadata format is:

```
{
'gaia' : {
         'parameters' : names of parameters (i.e. column headers) that exist in that
          survey
         'errors' : names of errors (i.e. column headers) for these parameters that
          exist in that survey
         'notes' : str, any notes to include on this parameter/error/etc.
         }
}

etc. for any supported survey
```

<br>

If a key is provided that is not the name of a supported survey, that key will be interpreted as a custom entry.

In this case, an additional 'values' key must be included, and the values/errors must be passed manually.

```
{
'custom' : {
           'parameters' : names of parameters
           'values' : parameter values
           'errors' : error values
           'notes' : str, any notes to include on this parameter/error/etc.
           }
}
```

**Returns**: dict

```
{
'type' : 'metadata'
'source' : int/str, source used to get data (None if a pos was used)
'pos': [ra,dec], pos used to get data (None if a source was used)
'plot' : bokeh figure object, the actual plot element
}
```

**NOTE:** Again, see the datapage_creation example as noted [above](#8-datapage-tools) for an example.

<br><br>

## 9. File Handling Tools<a id="9-file-handling-tools"></a>
These tools can be used to transform many returned ATK data structures into local files, and vice versa. This process is designed to be completely lossless, allowing for the exact data structure that was used to create a file to be recreated at a later date.

<br>

### 9.1. savedata<a id="91-savedata"></a>

This tool allows for the saving of many ATK data structures into local files.

Supported ATK data types: data, phot, bulkphot, image, sed, spectra, lightcurve

<br>

**Usage:**
```
savedata(data)
```

where:
```
data = any of the supported data structures from ATK
```

<br>

**Returns:** str

Name of file that the data was saved to

<br>

**Note:** Created files will contain a tag (e.g. ATKimage for an image). This is necessary, as it tells the readdata function (see below) how to interpret the file to recreate the original data structure.

<br>

### 9.2. readdata<a id="92-readdata"></a>

This tool allows for the lossless recreation of ATK data stuctures from files created using [savedata](#91-savedata).

<br>

**Usage:**

```
readdata(filename)
```

where:
```
filename = str, name/path of ATK file to read
```

<br>

**Returns:** ATK datastructure that was used to create the file

<br>

**Example:**
ATLAS lightcurve queries can take a long time (5-10 minutes). It is therefore useful to be able to save these to a local file for later use.

The following example shows how to fetch ATLAS lightcurve data, save it to a local file, read this local file to recover the data, and then plot it.
```
from AstroToolkit.Tools import lightcurvequery,savedata,readdata,plotlightcurve,showplot

lightcurve_data=lightcurvequery(source=6050296829033196032,survey='atlas')
fname=savedata(lightcurve_data)

recreated_lightcurve_data=readdata(fname)
showplot(plotlightcurve(recreated_lightcurve_data,colours=['red','orange','yellow']))
```

![](https://github.com/WD-planets/AstroToolkit/raw/latest/README_Images/struct_recreation_example.png?raw=true)

**Note:** The above example assumes that you have already set a default ATLAS username / password in the [config](#configuration).

<br>

### 9.3. export<a id="93-export"></a>

This allows for ATK figures to be saved as a .png file.

<br>

**Usage:**

```
export(plot,keephtml=True)
```

where:
```
plot = bokeh figure object
```
```
keephtml = bool. To get the filename for the png, the plot must first be saved as a regular bokeh static html file. This tells ATk whether to delete this file or keep it.
```

<br><br>

### 10. Miscellaneous Tools<a id="10-miscellaneous-tools"></a>

These are tools that do not fit into the above categories, but can still be useful for quickly performing some basic functions.

<br>

### 10.1. showplot<a id="101-showplot"></a>

This allows for the plot dictionaries returned by ATK to be opened in the browser. Doing so will also save the plot locally as a .html file. Can also be used with the same functionality as show() from bokeh for non-ATK bokeh figures.

<br>

**Usage:**
```
showplot(plot)
```

where:
```
plot = ATK plot dictionary
```

<br>

### 10.2. saveplot<a id="102-saveplot"></a>

This allows for the plot dictionaries returned by ATK to be saved to local .html files without opening them in the browser. Can also be used with the same functionality as save() from bokeh for non-ATK bokeh figures.

<br>

**Usage:**
```
saveplot(plot)
```

where:
```
plot = ATK plot dictionary
```

<br>

### 10.3. correctpm<a id="103-correctpm"></a>

Corrects for proper motion of an object between an input time and a target time.

<br>

**Usage:**

```
correctpm(inputtime,targettime,ra,dec,pmra,pmdec)
```

where:

```
inputtime = [year,month] where both entries are integers
```
```
targettime = [year,month] where both entries are integers
```
```
ra = float, ra of object in degrees
```
```
dec = float, dec of object in degrees
```
```
pmra = float, proper motion in ra direction of object in mas/yr
```
```
pmdec = float, proper motion in dec direction of object in mas/yr
```

<br>

**Returns:** list

[ra,dec] of object in degrees, corrected for proper motion.

<br>

**Example:**
To correct an object for proper motion to the year 2000:
```
from AstroToolkit.Tools import dataquery,correctpm

gaia_data=dataquery(survey='gaia',source=6050296829033196032)['data']
ra,dec,pmra,pmdec=gaia_data['ra'][0],gaia_data['dec'][0],gaia_data['pmra'][0],gaia_data['pmdec'][0]

ra_corrected,dec_corrected=correctpm([2016,0],[2000,0],ra,dec,pmra,pmdec)
```

<br>

### 10.4. getdistance<a id="104-getdistance"></a>

Simply calculates the distance (1/parallax) of an object given its parallax in mas (the unit of parallax in Gaia).

<br>

**Usage:**

```
getdistance(parallax)
```

where:
```
parallax = float, parallax of object in mas
```

<br>

**Returns:** float

distance of object in pc

<br>

### 10.5. convfromdeg<a id="105-convfromdeg"></a>

Converts deg coordinates to HMS/DMS format

<br>

**Usage:**

```
convfromdeg(pos)
```

where:
```
pos = list, [ra,dec] of object in degrees
```

<br>

**Returns:** list

[ra,dec] of object in [HMS,DMS]

<br>

### 10.6. convtodeg<a id="106-convtodeg"></a>

Converts HMS/DMS coordinates to deg.

<br>

**Usage:**

```
convfromdeg(pos)
```

where:
```
pos = list, [ra,dec] of object in [HMS/DMS] format, i.e. [[H,M,S],[D,M,S]]
```

<br>

**Returns:** list

[ra,dec] of object in deg

<br>

### 10.7. getcolumn<a id="107-getcolumn"></a>

Returns a list of values from a .fits file given a column name.

**Usage:**

```
getcolumn(file_name,col_name)
```

where:
```
file_name = str, name of file with or without the .fits extension.
```
```
col_name = str, name of column in file
```

**Returns:** list

list of values for that column

<br>

### 10.8. getsources<a id="108-getsources"></a>

Returns a list of source values from a .fits file given a column name, which can either be set in the tool or in the config.

**Usage:**

```
getsources(file_name,col_name=CONFIG)
```

where:
```
file_name = str, name of file with or without the .fits extension.
```
```
col_name = str, name of column containing Gaia sources
```

**Returns:** list

list of Gaia sources

<br>

### 10.9. getpositions<a id="109-getpositions"></a>

Returns a list of ra,dec values from a .fits file given ra and dec column names, which can either be set in the tool or in the config.

**Usage:**

```
getsources(file_name,col_name=CONFIG)
```

where:
```
file_name = str, name of file with or without the .fits extension.
```
```
col_name = str, comma separated string containing names of ra and dec columns. e.g. col_name = 'RA_ICRS,'DE_ICRS'. Must be in this order.
```

**Returns:** list

list of positions [[ra1,dec1],[ra2,dec2]], etc.