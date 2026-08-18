# SHIVA
SHIVA (a Simple and Helpful Interface for Variability Analysis) is a GUI wrapper for PARVATI (Profiles Analysis and Radial Velocities using Astronomical Tools for Investigation).
PARVATI may be either downloaded from [GitHub](https://github.com/mrainer74/parvati) or installed using pip:
```
python -m pip install parvati
```

> [!IMPORTANT]
> SHIVA >= v2 works only with PARVATI >= 2.0.0, previous versions are incompatible with the new `fit_profile` function of PARVATI.
> SHIVA >= v3 moved from Tkinter to PyQT6

## Introduction
SHIVA allows to use PARVATI with a detailed GUI where the user can set all the relevant parameters for the functions of PARVATI. It is possible to plot all the steps of the analysis, and everything is saved as FITS files with exhaustive headers.

SHIVA requires the following Python packages:
```
parvati >= 2.0.0
numpy >= 1.26.x
astropy >= 7.x.x
matplotlib >= 3.10
PyQt6
```

SHIVA has been tested only on Linux (Ubuntu 24.04) so far.
It is suggested to run SHIVA from a virtual environment. It does not need to be installed, it may be run with a simple command:
```
python shiva.py
```

SHIVA is organised in several tabs: a complete analysis from the reduced spectra to the profiles analysis should follow the order of the tabs. It is possible to work on a single spectrum or on all the spectra contained in a single directory.
All the input windows show tooltips with basic information on the input values.

## Left Tab 1: Normalisation
By default, SHIVA look for spectra in the same directory from where it is running, and it will save all the output in a new directory named `shiva_output`. Both the input and output directories may be specified in this tab. The `File/Pattern` entry allows to select a single file to work on or a pattern (using * as a wildcard), in which case SHIVA will work on all the files with the pattern in their names inside the input directory.
The spectra may be either ASCII or FITS files. The allowed formats are the same required by the `read_spectrum` function of PARVATI:
- monodimensional FITS files with the flux as the hdu[0].data and the wavelength in the hdu[0].header (CRVAL1, CDELT1, NAXIS1)
- FITS files in the e2ds format of HARPS/HARPS-N/SOPHIE, with the echelle ordes still unmerged and the wavelength information in the header in the *DRS CAL TH DEG LL and *DRS CAL TH COEFF LLXX keywords
- FITS tables with all the data in different fields of hdu[1].data OR in different hdus. By default, the wavelength will be read in the first field/hdu and the flux in the second field/hdu, but the number of the field/hdu may be specified. If there are any additional data as S/N and/or echelle order number and/or normalised flux and/or absolute errors, they may be specified here. If given, the S/N supersedes the errors, otherwise the errors will be transformed in S/N (S/N=flux/errors).
- ASCII files with at least two columns (wavelength and flux), but additional columns with S/N and/or echelle order number and/or normalised flux and/or absolute errors may be specified here. If given, the S/N supersedes the errors, otherwise the errors will be transformed in S/N (S/N=flux/errors). 
The normalised spectra are saved as FITS table with all the original information plus the normalised flux stored as fields in the hdu[1].data.

The numbers of columns/fields/hdus may be manually specified in the correspective entries OR the instrumentf may be selected from the scroll-down menu, and then the columns/fields/hdus are selected automatically.
If the instrument is `UNDEF` then the data will be searched with the default options: monodimensional FITS file with CRVAL1, CDELT1, NAXIS1, e2ds FITS file from HARPS/HARPS-N/SOPHIE, FITS table with wavelength and flux in the field/hdu 1 and 2, ASCII file with wavelength and flux in columns 1 and 2.

> [!NOTE] 
> - the numbers of the ASCII columns and the FITS fields start with 1, not 0
> - instead, the number of the hdu is the correct one: in the case of different hdu, the Primary hdu[0] is always empty
> - always specify the right wavelength unit: [a]ngstroms, [n]anometers or [m]icrons
> - when working with merged echelle spectra, use the `Subsets` option to achieve a good result
> - when working only with single line extraction, the profile normalisation may be enough, so this step may be skipped or used with the parameters `Subsets=0` and `Degree=0`.
> [!TIP]
> Read the GIANO-B ms1d data with the options: Wave=2, Flux=3, S/N=4, Orders=1
> Read the ESPRESSO S1D data with the options: Wave=1, Flux=3, Errors=4
> Read the ESPRESSO S2D data with the options: Wave=4, Flux=1, Errors=2
> Read the CARMENES (VIS and NIR) data with the options: wavecol=4, fluxcol=1, errcol=3

## Left Tab 2: Line Profile
This tab manages the extraction of single spectroscopic lines or the creation of mean line profiles using either the LSD or CCF methods. See the PARVATI README file for information on the functions `extract_line`, `compute_lsd` and `compute_ccf`.
The input and output directories and file/pattern value will be automatically updated from the normalisation parameters, but they may still be changed manually.
When creating a mean line profile, different kinds of masks may be used: ASCII files (VALD stellar mask, simple 2-columns file, normalised spectrum/model) or FITS files (mask or normalised spectrum/model. The CCF may also be computed weighting the contribution of the spectra using their S/N at the various wavelengths.
The resulting profiles are saved as FITS tables.

## Left Tab 3: Line Profile Analysis
The resulting profiles may be analysed using several PARVATI functions. SHIVA requires the profiles to be in the FITS format created in the [previous tab](#tab-2-line-profile).
The input/output folders and file pattern are used by all the functions in this tab.

### Left Tab 3a: Profile Normalisation
The profiles may be normalised with a simple linear fitting of the continuum, defined as the region outside the `RV min` and `RV max` parameters. 
If more than one profile is given as input, then the `St. Dev.` option results in the computation of an average line profile and the standard deviation of all the profiles from the mean. They will saved in a `line_mean_std.txt` file.
The resulting profiles are saved as FITS tables.
> [!TIP]
> Using the `St. Dev.` option allows a quick look at the impact and location of any line profile variation.
>[!NOTE]
>The output suffix of this step will automatically update the input file pattern, as the subsequent analysis should be carried out on the normalised profiles.

### Left Tab 3b: Profile Fitting
It is possible to fit up to 3 independent components (multiple spectroscopic system), choosing between all the available functions in PARVATI: Gaussian, Asymmetric Gaussian, Supergaussian, Lorentzian, Voigt, Rotational profile. If the guess RV value is outside the RV range of the profile, it will be automatically shifted to the position of the minimum flux of the profile when running the fit.
Using the absolute errors of the flux to perform the fit will result in smaller errors on the fit if the flux errors are reliable, otherwise it is suggested to uncheck the `Use errors` option.
The resulting fit profile is saved as an additional column in the input FITS profile, along with the fitting parameters values that are also saved by updating the header content.
The fit report may be saved as an ASCII file by checking the `Save fit report` option.
> [!TIP]
> Always use either the Gaussian or rotational fit if you plan to perform also the subsequent analysis steps (moments, bisector and Fourier Transform). These fits will allow to better define the line limits.

### Left Tab 3c: Line Analysis: Moments/Bisector/Fourier Transform
The same line limits may be set for the subsequent analysis grouped in the same tab.
The line limits must be defined, either by inputting two fixed RV values (lower and upper limit) or by selecting the `Gaussian` or `Rotational` parameter: when the `Gaussian` limits are chosen, the line limits are defined as the Gaussian RV values +/- 3 sigma, while when the `Rotational` limits are chosen, the line limits are defined as the rotational RV values +/- *v*sin*i*.

#### Moments
The first 5 line moments are computed (from m0 to m4), and the skewness and the kurtosis are derived from m3 and m4. 
The resulting values are saved in the input profile FITS files, updating only the header content.

#### Bisector
The line bisector and the bisector's span are computed.
The resulting values are saved in the input profile FITS files, updating only the header content.

#### Fourier Transform
The Fourier Transform (FT) of the symmetrised line is computed. The symmetrisation process yields another RV estimation, while the positions of the first 3 zeroes of the FT results in 3 estimation of the stellar *v*sin*i*, and average value and an indicator of differential rotation.
The resulting FTs are saved as FITS tables, and the relevant output values are also saved in the input profile FITS files, updating only the header content.

## Left Tab 4: Time Series
Once the previous analysis have been performed, it is possible to plot the resulting time series and save them as text files.

## Rght Tabs: Plot and Log
The right side of SHIVA display the plots generated by all the tabs and a log of the process. The log may be saved as a text file using the `Save log` button.

## Test files
A few reduced echelle spectra and two VALD stellar masks are given in the `tests` directory, to help familiarising with SHIVA.

