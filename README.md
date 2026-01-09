# SHIVA
SHIVA (a Simple and Helpful Interface for Variability Analysis) is a GUI wrapper for PARVATI (Profiles Analysis and Radial Velocities using Astronomical Tools for Investigation).
PARVATI may be either downloaded from [GitHub](https://github.com/mrainer74/parvati) or installed using pip:
```
python -m pip install parvati
```

## Introduction
SHIVA allows to use PARVATI with a detailed GUI where the user can set all the relevant parameters for the functions of PARVATI. It is possible to plot all the steps of the analysis, and everything is saved as FITS files with exhaustive headers.

SHIVA requires the following Python packages:
```
parvati >= 1.0.5
numpy >= 1.26.x
astropy >= 7.x.x
matplotlib >= 3.10
threading
python-tk : same version as python
TkToolTip >= 1.2
```

SHIVA has been tested on Linux (Ubuntu 24.04) and MacOS 26 (Tahoe) so far.
It is suggested to run SHIVA from a virtual environment. It does not need to be installed, it may be run with a simple command:
```
python shiva.py
```

SHIVA is organised in several tabs: a complete analysis from the reduced spectra to the profiles analysis should follow the numbering of the tabs. It is possible to work on a single spectrum or on all the spectra contained in a single directory.
All the input windows show tooltips with basic information of the input values.

## Tab 1: Normalisation
By default, SHIVA look for spectra in the same directory from where it is running, and it will save all the output in a new directory named `shiva_output`. Both the input and output directories may be specified in this tab. The `File/Pattern` entry allows to select a single file to work on or a pattern (using * as a wildcard), in which case SHIVA will work on all the files with the pattern in their names inside the input directory.
The spectra may be either ASCII or FITS files. The allowed formats are the same required by the `read_spectrum` function of PARVATI:
- monodimensional FITS files with the flux as the hdu[0].data and the wavelength in the hdu[0].header (CRVAL1, CDELT1, NAXIS1)
- FITS tables with all the data in hdu[1].data. By default, the wavelength will be read in the first field and the flux in the second field, but the number of the field may be specified. If there are any additional data as S/N and/or echelle order number and/or normalised flux and/or absolute errors, they may be specified here. If given, the S/N supersedes the errors, otherwise the errors will be transformed in S/N (S/N=flux/errors).
- ASCII files with at least two columns (wavelength and flux), but additional columns with S/N and/or echelle order number and/or normalised flux and/or absolute errors may be specified here. If given, the S/N supersedes the errors, otherwise the errors will be transformed in S/N (S/N=flux/errors). 
The normalised spectra are saved as FITS table with all the original information plus the normalised flux stored as fields in the hdu[1].data.

[!NOTE] 
- the numbers of the ASCII columns and the FITS fields start with 1, not 0
- always specify the right wavelength unit: [a]ngstroms, [n]anometers or [m]icrons
- when working with merged echelle spectra, use the `Subsets` option to achieve a good result
- when working only with single line extraction, the profile normalisation may be enough, so this step may be skipped or used with the parameters `Subsets=0` and `Degree=0`.

## Tab 2: Line Profile
This tab manages the extraction of single spectroscopic lines or the creation of mean line profiles using either the LSD or CCF methods. See the PARVATI README file for information on the functions `extract_line`, `compute_lsd` and `compute_ccf`.
The input and output directories and file/pattern value will be automatically updated from the normalisation parameters, but they may still be changed manually.
When creating a mean line profile, different kinds of masks may be used: ASCII files (VALD stellar mask, simple 2-columns file, normalised spectrum/model) or FITS files (mask or normalised spectrum/model. The CCF may also be computed weighting the contribution of the spectra using their S/N at the various wavelengths.
The resulting profiles are saved as FITS tables.

## Tab 3: Line Analysis
The resulting profiles may be analysed using several PARVATI functions. SHIVA requires the profiles to be in the FITS format created in the [previous tab](#tab-2-line-profile).

### Tab 3a: Profile Normalisation
The profiles may be normalised with a simple linear fitting of the continuum, defined as the region outside the `RV min` and `RV max` parameters. 
If more than one profile is given as input, then the `St. Dev.` option results in the computation of an average line profile and the standard deviation of all the profiles from the mean. This will saved in a `line_mean_std.txt` file (data), and a `line_mean_std.png` file (plot).
The resulting profiles are saved as FITS tables.
[!TIP]
Using the `St. Dev.` option allows a quick look at the impact and location of any line profile variation.

### Tab 3b: Profile Fitting
By default, all 4 possible functions (Gaussian, Lorentzian, Voigt, rotational) are used. If the guess RV value is outside the RV range of the profile, it will be automatically shifted to the middle of the profile when running the fit.
Using the absolute errors of the flux to perform the fit will result in smaller errors on the fit if the flux errors are reliable, otherwise it is suggested to uncheck the `Use errors` option.
The resulting fit profiles are saved as FITS tables, and the fitting parameters values are also saved in the input profile FITS files, updating only the header content.
[!TIP]
Always use either the Gaussian or rotational fit if the subsequent analysis steps (moments, bisector and Fourier Transform) are done. These fits will allow to better define the line limits.

### Tab 3c: Moments
The first 5 line moments are computed (from m0 to m4), and the skewness and the kurtosis are derived from m3 and m4. 
The line limits must be defined, either by inputting two fixed RV values (lower and upper limit) or by checking the `Gaussian` or `Rotational` box: when the `Gaussian` box is checked the line limits are defined as the Gaussian RV values +/- 3 sigma, while when the `Rotational` box is checked the line limits are defined as the rotational RV values +/- *v*sin*i*.
The resulting values are saved as FITS tables, and the moments values are also saved in the input profile FITS files, updating only the header content.

### Tab 3d: Bisector
The line bisector and the bisector's span are computed.
The line limits must be defined, either by inputting two fixed RV values (lower and upper limit) or by checking the `Gaussian` or `Rotational` box: when the `Gaussian` box is checked the line limits are defined as the Gaussian RV values +/- 3 sigma, while when the `Rotational` box is checked the line limits are defined as the rotational RV values +/- *v*sin*i*.
The resulting values are saved as FITS tables, and the moments values are also saved in the input profile FITS files, updating only the header content.

### Tab 3e: Fourier Transform
The Fourier Transform (FT) of the symmetrised line is computed. The symmetrisation process yields another RV estimation, while the positions of the first 3 zeroes of the FT results in 3 estimation of the stellar *v*sin*i*, and average value and an indicator of differential rotation.
The line limits must be defined, either by inputting two fixed RV values (lower and upper limit) or by checking the `Gaussian` or `Rotational` box: when the `Gaussian` box is checked the line limits are defined as the Gaussian RV values +/- 3 sigma, while when the `Rotational` box is checked the line limits are defined as the rotational RV values +/- *v*sin*i*.
The resulting FTs are saved as FITS tables, and the relevant output values are also saved in the input profile FITS files, updating only the header content.

### Other tabs: plot and log
The middle tabs of SHIVA display the plots generated by all the tabs (if the `Plot` parameter is checked) and a log of the process. The log may be saved as a text file using the `Save log` button.

## Test files
A few reduced echelle spectra and two VALD stellar masks are given in the `tests` directory, to help familiarising with SHIVA.

