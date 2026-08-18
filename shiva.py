"""
SHIVA PyQt GUI wrapper for PARVATI.
This module provides a PyQt6-based user interface for the SHIVA workflow.

SHIVA: a Simple and Helpful Interface for Variability Analysis
SHIVA is a PyQt GUI wrapper for PARVATI: Profiles Analysis and Radial 
Velocities using Astronomical Tools for Investigation 
(a Python package to compute and analyse stellar mean line profiles)

Written by Monica Rainer
with the help of VSCode AI to convert from the SHIVA v2 Tkinter GUI to PyQt

    SHIVA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY. 
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""

import os, sys
import datetime
import glob
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', AstropyWarning)

import matplotlib
matplotlib.use('QtAgg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT

import parvati as pa

from PyQt6 import QtCore, QtWidgets
from PyQt6.QtCore import Qt

__version__ = '3.1.0'

class Var:
    def __init__(self, value=None):
        self._value = value
        self._callbacks = []

    def set(self, value):
        self._value = value
        for callback in self._callbacks:
            callback(value)

    def get(self):
        return self._value

    def trace(self, mode, callback):
        self.trace_add(mode, callback)

    def trace_add(self, mode, callback):
        if callback not in self._callbacks:
            self._callbacks.append(callback)

    def __str__(self):
        return str(self._value)

    def __repr__(self):
        return f"Var({self._value!r})"

class Worker(QtCore.QRunnable):
    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def run(self):
        self.fn()

class ShivaQtApp(QtWidgets.QMainWindow):
    log_signal = QtCore.pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.basedir = os.getcwd()
        self.version = ''.join(('SHIVA v', __version__))
        try:
            self.pa_version = ''.join(('PARVATI v', pa.__version__))
        except AttributeError:
            self.pa_version = 'PARVATI vXXX'
        self.key_prg = 'SP'
        self.sp_logo = 'sp_logo.png'
        self.outdir = 'shiva_output'
        self.back = '#d5e6f5' # very light blue
        self.color = '#cf0202' # dark red
        self.runbutton = '#9bfaa1' # very light green
        self.browse = '#f9fac3' # very light yellow 
        self.define_keys()
        self.define_entries()
        self.set_entries()
        self.define_plot()
        self.setWindowTitle('SHIVA - Simple and Helpful Interface for Variability Analysis using PARVATI')
        self.resize(1600, 900)
        self.init_ui()
        self.log_signal.connect(self.append_log)

    def define_keys(self):
        # General
        self.input_jd = ['MJD-OBS']
        self.key_ver = (f'HIERARCH {self.key_prg} VERSION', 'SHIVA version used')
        self.pa_ver = (f'HIERARCH PARVATI VERSION', 'PARVATI version used')
        
        # Normalisation
        self.key_nor = (f'HIERARCH {self.key_prg} NOR DEG', 'Polynomial degree for normalisation.')
        self.key_sub = (f'HIERARCH {self.key_prg} NOR SUBSETS', 'Subsets independently normalised.')
        self.key_refine = (f'HIERARCH {self.key_prg} NOR REFINE', 'Refined normalisation')
        self.key_jd = (f'HIERARCH {self.key_prg} JD', 'JD value in input spectrum ')
        self.key_swave = (f'HIERARCH {self.key_prg} WAVE UNITS', 'Wavelength unit (a=Angstrom, n=nanometer, m=micron)')
        
        # Line Profile
        self.base_key = {'lsd': f'HIERARCH {self.key_prg} LSD', 'ccf': f'HIERARCH {self.key_prg} CCF', 'line': f'HIERARCH {self.key_prg} EXT'}
        self.key_lsdin = ('INPUT', 'Input spectrum used to compute profile')
        self.key_mask = ('MASK', 'Mask used to compute profile')
        self.key_mask_invert = ('MASK INVERT', 'Absorption mask, inverted before being used')
        self.key_mask_spectrum = ('MASK SPECTRUM', 'Spectrum (observed or template) used as mask - with continuum')
        self.key_els = ('MASK ELS', 'Mask elements used')
        self.key_noels = ('MASK NO ELS', 'Mask elements excluded')
        self.key_dlow = ('MASK DLOW', 'Mask minimum depth used')
        self.key_dup = ('MASK DUP', 'Mask maximum depth used')
        self.key_wmin = ('MASK WMIN', 'Mask minimum wavelength used')
        self.key_wmax = ('MASK WMAX', 'Mask maximum wavelength used')
        self.key_balmer = ('MASK BALMER', 'Balmer regions excluded')
        self.key_tell = ('MASK TELLURIC', 'Telluric regions excluded')
        self.key_cosmic = ('COSMIC', 'Cosmic rays removed before the profile')
        self.key_clean = ('CLEAN', 'Clean profile (smoothing applied)')
        self.key_rvmin = ('RVMIN', 'Minimum RV of the profile (km/s)')
        self.key_rvmax = ('RVMAX', 'Maximum RV of the profile (km/s)')
        self.key_rvstep = ('RVSTEP', 'RV step of the profile (km/s)')
        self.key_ccfweight = (' '.join((self.base_key['ccf'], 'WEIGHTED')), 'Weighted CCF, using the S/N')
        
        # Line Analysis
        self.key_lainp = (f'HIERARCH {self.key_prg} LA INPUT', 'Input file for the line analysis')
        self.key_lnor = (f'HIERARCH {self.key_prg} PRF NOR', 'Normalised profile')
        self.key_lnrvmin = (f'HIERARCH {self.key_prg} PRF NOR RVMIN', 'Minimum RV for continuum definition (km/s)')
        self.key_lnrvmax = (f'HIERARCH {self.key_prg} PRF NOR RVMAX', 'Maximum RV for continuum definition (km/s)')
        
        # Line fit
        self.key_comps = (f'HIERARCH {self.key_prg} FIT COMP', 'Number of components fitted')
        self.key_res = (f'HIERARCH {self.key_prg} FIT RES', 'Instrumental resolution used in the fit')
        
        self.key_fit1 = (f'HIERARCH {self.key_prg} FIT1', 'Fit function for component 1')
        self.key_fit2 = (f'HIERARCH {self.key_prg} FIT2', 'Fit function for component 2')
        self.key_fit3 = (f'HIERARCH {self.key_prg} FIT3', 'Fit function for component 3')
        
        self.key_rvguess1 = (f'HIERARCH {self.key_prg} FIT1 RV GUESS', 'RV (km/s) guess value for component 1')
        self.key_rvguess2 = (f'HIERARCH {self.key_prg} FIT2 RV GUESS', 'RV (km/s) guess value for component 2')
        self.key_rvguess3 = (f'HIERARCH {self.key_prg} FIT3 RV GUESS', 'RV (km/s) guess value for component 3')
        
        self.key_wguess1 = (f'HIERARCH {self.key_prg} FIT1 WIDTH GUESS', 'Line width (km/s) guess value for component 1')
        self.key_wguess2 = (f'HIERARCH {self.key_prg} FIT2 WIDTH GUESS', 'Line width (km/s) guess value for component 2')        
        self.key_wguess3 = (f'HIERARCH {self.key_prg} FIT3 WIDTH GUESS', 'Line width (km/s) guess value for component 3')   
             
        self.key_fitld1 = (f'HIERARCH {self.key_prg} FIT1 LD', 'Linear Limb Darkening for component 1')
        self.key_fitld2 = (f'HIERARCH {self.key_prg} FIT2 LD', 'Linear Limb Darkening for component 2')
        self.key_fitld3 = (f'HIERARCH {self.key_prg} FIT3 LD', 'Linear Limb Darkening for component 3')

        self.key_rv1 = (f'HIERARCH {self.key_prg} FIT1 RV', 'RV (km/s) fit result for component 1')
        self.key_rv2 = (f'HIERARCH {self.key_prg} FIT2 RV', 'RV (km/s) fit result for component 2')
        self.key_rv3 = (f'HIERARCH {self.key_prg} FIT3 RV', 'RV (km/s) fit result for component 3')
        
        self.key_rverr1 = (f'HIERARCH {self.key_prg} FIT1 RVERR', 'RV (km/s) fit error for component 1')
        self.key_rverr2 = (f'HIERARCH {self.key_prg} FIT2 RVERR', 'RV (km/s) fit error for component 2')
        self.key_rverr3 = (f'HIERARCH {self.key_prg} FIT3 RVERR', 'RV (km/s) fit error for component 3')
        
        self.key_fwhm1 = (f'HIERARCH {self.key_prg} FIT1 WIDTH', 'FWHM/vsini (km/s) fit for component 1')
        self.key_fwhm2 = (f'HIERARCH {self.key_prg} FIT2 WIDTH', 'FWHM/vsini (km/s) fit for component 2')
        self.key_fwhm3 = (f'HIERARCH {self.key_prg} FIT3 WIDTH', 'FWHM/vsini (km/s) fit for component 3')
        
        self.key_fwhmerr1 = (f'HIERARCH {self.key_prg} FIT1 FWHMERR', 'FWHM/vsini (km/s) fit error for component 1')
        self.key_fwhmerr2 = (f'HIERARCH {self.key_prg} FIT2 FWHMERR', 'FWHM/vsini (km/s) fit error for component 2')
        self.key_fwhmerr3 = (f'HIERARCH {self.key_prg} FIT3 FWHMERR', 'FWHM/vsini (km/s) fit error for component 3')
        
        self.key_fitew1 = (f'HIERARCH {self.key_prg} FIT1 EW', 'EW (km/s) fit for component 1')
        self.key_fitew2 = (f'HIERARCH {self.key_prg} FIT2 EW', 'EW (km/s) fit for component 2')
        self.key_fitew3 = (f'HIERARCH {self.key_prg} FIT3 EW', 'EW (km/s) fit for component 3')
        
        self.key_fitewerr1 = (f'HIERARCH {self.key_prg} FIT1 EWERR', 'EW (km/s) fit error for component 1')
        self.key_fitewerr2 = (f'HIERARCH {self.key_prg} FIT2 EWERR', 'EW (km/s) fit error for component 2')
        self.key_fitewerr3 = (f'HIERARCH {self.key_prg} FIT3 EWERR', 'EW (km/s) fit error for component 3')
        
        # Moments
        self.key_momlim = (f'HIERARCH {self.key_prg} MOM LIMITS', 'Line limits for moments computation')
        self.key_mom0 = (f'HIERARCH {self.key_prg} MOM M0', '0th moment (EW)')
        self.key_mom0err = (f'HIERARCH {self.key_prg} MOM M0ERR', 'Error on 0th moment')
        self.key_mom1 = (f'HIERARCH {self.key_prg} MOM M1', '1st moment (RV)')
        self.key_mom1err = (f'HIERARCH {self.key_prg} MOM M1ERR', 'Error on 1st moment')
        self.key_mom2 = (f'HIERARCH {self.key_prg} MOM M2', '2nd moment (Variance)')
        self.key_mom2err = (f'HIERARCH {self.key_prg} MOM M2ERR', 'Error on 2nd moment')
        self.key_mom3 = (f'HIERARCH {self.key_prg} MOM M3', '3rd moment (seed for skewness)')
        self.key_mom3err = (f'HIERARCH {self.key_prg} MOM M3ERR', 'Error on 3rd moment')
        self.key_mom4 = (f'HIERARCH {self.key_prg} MOM M4', '4th moment (seed for kurtosis)')
        self.key_mom4err = (f'HIERARCH {self.key_prg} MOM M4ERR', 'Error on 4th moment')
        self.key_momfwhm = (f'HIERARCH {self.key_prg} MOM FWHM', 'FWHM from 2nd moment')
        self.key_momfwhmerr = (f'HIERARCH {self.key_prg} MOM FWHMERR', 'Error on FWHM from 2nd moment')
        self.key_momskew = (f'HIERARCH {self.key_prg} MOM SKEWNESS', 'Skewness from 3rd moment')
        self.key_momskewerr = (f'HIERARCH {self.key_prg} MOM SKERR', 'Error on Skewness from 3rd moment')
        self.key_momkurt = (f'HIERARCH {self.key_prg} MOM KURTOSIS', 'Kurtosis from 4th moment')
        self.key_momkurterr = (f'HIERARCH {self.key_prg} MOM KURERR', 'Error on Kurtosis from 4th moment')
        
        # Bisector
        self.key_bislim = (f'HIERARCH {self.key_prg} BIS LIMITS', 'Line limits for bisector computation (km/s)')
        self.key_bispan = (f'HIERARCH {self.key_prg} BIS SPAN', 'Bisector span (km/s)')
        self.key_biserr = (f'HIERARCH {self.key_prg} BIS SPANERR', 'Bisector span error (km/s)')
        
        # Fourier
        self.key_foulim = (f'HIERARCH {self.key_prg} FOU LIMITS', 'Line limits for Fourier Transform (km/s)')
        self.key_fouerr = (f'HIERARCH {self.key_prg} FOU ERROR', 'Mean error of the FTT power')
        self.key_fouz1 = (f'HIERARCH {self.key_prg} FOU Z1', 'First zero position')
        self.key_fouz1err = (f'HIERARCH {self.key_prg} FOU Z1ERR', 'Error on first zero position')
        self.key_fouz2 = (f'HIERARCH {self.key_prg} FOU Z2', 'Second zero position')
        self.key_fouz2err = (f'HIERARCH {self.key_prg} FOU Z2ERR', 'Error on second zero position')
        self.key_fouz3 = (f'HIERARCH {self.key_prg} FOU Z3', 'Third zero position')
        self.key_fouz3err = (f'HIERARCH {self.key_prg} FOU Z3ERR', 'Error on third zero position')
        self.key_fouratio = (f'HIERARCH {self.key_prg} FOU RATIO', 'Second/First zero ratio')
        self.key_fouratioerr = (f'HIERARCH {self.key_prg} FOU RATIOERR', 'Error on Second/First zero ratio')
        self.key_fouvsini1 = (f'HIERARCH {self.key_prg} FOU VSINI1', 'Vsini from first zero position')
        self.key_fouvsini1err = (f'HIERARCH {self.key_prg} FOU VSINI1ERR', 'Error on vsini from first zero position')
        self.key_fouvsini2 = (f'HIERARCH {self.key_prg} FOU VSINI2', 'Vsini from second zero position')
        self.key_fouvsini2err = (f'HIERARCH {self.key_prg} FOU VSINI2ERR', 'Error on vsini from second zero position')
        self.key_fouvsini3 = (f'HIERARCH {self.key_prg} FOU VSINI3', 'Vsini from third zero position')
        self.key_fouvsini3err = (f'HIERARCH {self.key_prg} FOU VSINI3ERR', 'Error on vsini from third zero position')
        self.key_fouvsini = (f'HIERARCH {self.key_prg} FOU VSINI', 'Mean vsini from all zero positions')
        self.key_fouvsinierr = (f'HIERARCH {self.key_prg} FOU VSINIERR', 'Error on mean vsini from first all positions')
        self.key_fourv = (f'HIERARCH {self.key_prg} FOU RV', 'RV derived from symmetrising the line prior to Fourier Tranform')

    def define_entries(self):
        self.nor_indir = Var()
        self.nor_spec = Var()
        self.option_instr = Var()
        self.units_default = Var()
        self.wave_frame = Var()
        self.wavecol = Var()
        self.fluxcol = Var()
        self.snrcol = Var()
        self.errcol = Var()
        self.echcol = Var()
        self.degree = Var()
        self.n_ord = Var()
        self.refine = Var()
        self.spec_unit = Var()
        self.spec_vacuum = Var()
        self.nor_outdir = Var()
        self.nor_output = Var()
        self.prf_indir = Var()
        self.prf_spec = Var()
        self.prf_wavecol = Var()
        self.prf_fluxcol = Var()
        self.prf_nfluxcol = Var()
        self.prf_snrcol = Var()
        self.prf_spec_unit = Var()
        self.prf_outdir = Var()
        self.ext_wave = Var()
        self.rvmin = Var()
        self.rvmax = Var()
        self.rvstep = Var()
        self.ext_output = Var()
        self.mask = Var()
        self.mask_invert = Var()
        self.mask_spectrum = Var()
        self.mask_cwave = Var()
        self.mask_cflux = Var()
        self.mask_unit = Var()
        self.mask_units_default = Var()
        self.mask_wave_frame = Var()
        self.mask_vacuum = Var()
        self.mask_dlow = Var()
        self.mask_dup = Var()
        self.mask_wmin = Var()
        self.mask_wmax = Var()
        self.mask_els = Var()
        self.mask_noels = Var()
        self.mask_balmer = Var()
        self.mask_tell = Var()
        self.cosmic = Var()
        self.clean = Var()
        self.ccfweight = Var()
        self.do_lsd = Var()
        self.do_ccf = Var()
        self.prf_output = Var()
        self.messages = Var()
        self.la_indir = Var()
        self.la_spec = Var()
        self.la_outdir = Var()
        self.option_fit = Var()
        self.option_fit2 = Var()
        self.option_fit3 = Var()
        self.option_limits = Var()
        self.limitlow = Var()
        self.limitup = Var()
        self.std = Var()
        self.norprf_output = Var()
        self.chosen_fit = Var()
        self.chosen_fit2 = Var()
        self.chosen_fit3 = Var()
        self.rv0 = Var()
        self.rv02 = Var()
        self.rv03 = Var()
        self.width = Var()
        self.width2 = Var()
        self.width3 = Var()
        self.ld = Var()
        self.ld2 = Var()
        self.ld3 = Var()
        self.resolution = Var()
        self.fit_errs = Var()
        self.save_report = Var()
        self.la_limitlow = Var()
        self.la_limitup = Var()
        self.fou_output = Var()        
        self.ts_indir = Var()
        self.ts_spec = Var()
        self.ts_outdir = Var()
        self.ts_comp = Var()
        self.ts_xplot = Var()
        self.ts_yplot = Var()
        self.ts_error = Var()
        self.ts_save = Var()
        self.thread_running = Var()
        self.abort_value = Var()

    def set_entries(self):
        self.option_values = ['DEFAULT', 'ESPRESSO S1D', 'ESPRESSO S2D', 'GIANO-B MS1D', 'CARMENES', 'HARPS(N) S1D NEW DRS', 'HARPS(N) S2D NEW DRS']
        self.units_values = ['angstroms', 'nm', 'micron']
        self.wave_values = ['Vacuum', 'Air']
        self.mask_units_values = ['angstroms', 'nm', 'micron']
        self.mask_wave_values = ['Vacuum', 'Air']
        self.option_la_fit = ['Gaussian', 'Rotational', 'Lorentzian', 'Voigt', 'Asymmetric Gaussian', 'Supergaussian']
        self.option_la_fits = ['None','Gaussian', 'Rotational', 'Lorentzian', 'Voigt', 'Asymmetric Gaussian', 'Supergaussian']
        self.option_la_values = ['Gaussian', 'Rotational', 'Manual']
        self.option_ts_plot = ['JD', 'RV', 'EW', 'width', 'm0_EW', 'm1_RV', 'm2_sigma', 'skewness', 'kurtosis', 'bispan', 'vsini_Fourier', 'q2/q1_Fourier']
        self.option_ts_comp = ['1', '2', '3']

        self.nor_indir.set(self.basedir)
        self.nor_spec.set('*.fits')
        self.option_instr.set(self.option_values[0])
        self.units_default.set(self.units_values[0])
        self.wave_frame.set(self.wave_values[0])        
        self.wavecol.set(1)
        self.fluxcol.set(2)
        self.snrcol.set(0)
        self.errcol.set(0)
        self.echcol.set(0)
        self.degree.set(2)
        self.n_ord.set(0)
        self.refine.set(0)
        self.spec_unit.set('a')
        self.spec_vacuum.set(1)
        self.nor_outdir.set(os.path.join(self.basedir, self.outdir))
        self.nor_output.set('_nor.fits')
        self.prf_indir.set(self.nor_outdir.get())
        self.prf_spec.set('*_nor.fits')
        self.prf_wavecol.set(1)
        self.prf_fluxcol.set(2)
        self.prf_nfluxcol.set(3)
        self.prf_snrcol.set(4)
        self.prf_spec_unit.set('a')
        self.prf_outdir.set(self.nor_outdir.get())
        self.ext_wave.set(6562.801)
        self.rvmin.set(-100)
        self.rvmax.set(100)
        self.rvstep.set(1)
        self.ext_output.set('_ext.fits')
        self.mask.set('')
        self.mask_invert.set(0)
        self.mask_spectrum.set(0)
        self.mask_cwave.set(1)
        self.mask_cflux.set(2)
        self.mask_unit.set('a')
        self.mask_units_default.set(self.mask_units_values[0])
        self.mask_wave_frame.set(self.mask_wave_values[0])
        self.mask_vacuum.set(1)
        self.mask_dlow.set(0.01)
        self.mask_dup.set(1)
        self.mask_wmin.set(0)
        self.mask_wmax.set(0)
        self.mask_els.set('')
        self.mask_noels.set('')
        self.mask_balmer.set(1)
        self.mask_tell.set(1)
        self.cosmic.set(0)
        self.clean.set(0)
        self.ccfweight.set(0)
        self.do_lsd.set(0)
        self.do_ccf.set(1)
        self.prf_output.set('_prf.fits')
        self.messages.set(f"{self.version}\n{self.pa_version}\n")
        self.la_indir.set(self.nor_outdir.get())
        self.la_spec.set('*_prf.fits')
        self.la_outdir.set(self.nor_outdir.get())
        self.option_fit.set(self.option_la_fit[0])
        self.option_fit2.set(self.option_la_fits[0])
        self.option_fit3.set(self.option_la_fits[0])
        self.chosen_fit.set(self.option_la_fit[0])
        self.chosen_fit2.set(self.option_la_fits[0])
        self.chosen_fit3.set(self.option_la_fits[0])
        self.option_limits.set(self.option_la_values[0])
        self.limitlow.set(-80)
        self.limitup.set(80)
        self.std.set(1)
        self.norprf_output.set('_pfn.fits')
        self.rv0.set(0)
        self.rv02.set(0)
        self.rv03.set(0)
        self.width.set(10)
        self.width2.set(10)
        self.width3.set(10)
        self.ld.set(0.6)
        self.ld2.set(0.6)
        self.ld3.set(0.6)
        self.resolution.set(0)
        self.fit_errs.set(1)
        self.save_report.set(0)
        #self.la_limitlow.set(-80)
        #self.la_limitup.set(80)
        self.fou_output.set('_fou.fits')
        self.ts_indir.set(self.nor_outdir.get())
        self.ts_spec.set('*_pfn.fits')
        self.ts_outdir.set(self.nor_outdir.get())
        self.ts_comp.set(self.option_ts_comp[0])
        self.ts_xplot.set(self.option_ts_plot[0])
        self.ts_yplot.set(self.option_ts_plot[1])
        self.ts_error.set(0)
        self.ts_save.set(0)
        self.thread_running.set(0)
        self.abort_value.set(False)

    def define_plot(self):
        self.do_plot = False
        self.xvalues = None
        self.yvalues = None
        self.y_add = None
        self.y_res = None
        self.ymin = None
        self.ymax = None
        self.limits = (None, None)
        self.logscale = False

    def init_ui(self):
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QtWidgets.QGridLayout()
        central_widget.setLayout(main_layout)

        
        self.tab_widget = QtWidgets.QTabWidget()
        self.tab_widget.setStyleSheet(f'QTabBar {{font: bold;}}')
        self.tab_widget.addTab(self.create_normalisation_tab(), 'Normalisation')
        self.tab_widget.addTab(self.create_profile_tab(), 'Line Profile')
        self.tab_widget.addTab(self.create_analysis_tab(), 'Line Profile Analysis')
        self.tab_widget.addTab(self.create_timeseries_tab(), 'Time Series')

        #left_panel = QVBoxLayout()   
        #left_panel.addWidget(self.create_log_group())
        #left_panel.setStretch(0, 1)

        right_panel = QtWidgets.QVBoxLayout()
        right_panel.addWidget(self.create_plot_group())
        right_panel.addWidget(self.create_log_group())
        #right_panel.addWidget(self.create_quit_group())
        right_panel.setStretch(0, 3)
        right_panel.setStretch(1, 2)

        main_layout.addWidget(self.tab_widget, 0, 0)
        #main_layout.addLayout(left_panel, 1, 0)
        main_layout.addLayout(right_panel, 0, 1)
        main_layout.setColumnStretch(0, 1)
        main_layout.setColumnStretch(1, 2)

    def create_normalisation_tab(self):
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
                
        # Input/output parameters
        inout_group = QtWidgets.QGroupBox('Input/Output')
        
        # input grid
        in_layout = QtWidgets.QGridLayout()
        
        input_norm_label = QtWidgets.QLabel("Input folder:")        
        self.nor_indir_edit = QtWidgets.QLineEdit(self.nor_indir.get())
        self.nor_indir_edit.setToolTip('Folder with input FITS/ASCII spectra')
        browse_dir = QtWidgets.QPushButton('Browse')
        browse_dir.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')
        browse_dir.clicked.connect(self.nor_load_indir)
        
        in_layout.addWidget(input_norm_label,0,0)
        in_layout.addWidget(self.nor_indir_edit,0,1)
        in_layout.addWidget(browse_dir,0,2)

        # File/pattern folder
        spec_norm_label = QtWidgets.QLabel("File/Pattern:")        
        self.nor_spec_edit = QtWidgets.QLineEdit(self.nor_spec.get())
        self.nor_spec_edit.setToolTip('Select pattern OR single FITS/ASCII spectrum')
        browse_file = QtWidgets.QPushButton('Browse')
        browse_file.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')
        browse_file.clicked.connect(self.nor_load_file)
        
        in_layout.addWidget(spec_norm_label,1,0)
        in_layout.addWidget(self.nor_spec_edit,1,1)
        in_layout.addWidget(browse_file,1,2)

        # Output folder
        out_norm_label = QtWidgets.QLabel("Output folder:")
        self.nor_outdir_edit = QtWidgets.QLineEdit(self.nor_outdir.get())
        self.nor_outdir_edit.setToolTip('Output folder (it will be created if needed)')        
        browse_out = QtWidgets.QPushButton('Browse')
        browse_out.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')
        browse_out.clicked.connect(self.nor_load_outdir)

        in_layout.addWidget(out_norm_label,2,0)
        in_layout.addWidget(self.nor_outdir_edit,2,1)
        in_layout.addWidget(browse_out,2,2)
        
        # Suffix
        sfx_norm_label = QtWidgets.QLabel("Output suffix:")
        self.nor_output_edit = QtWidgets.QLineEdit(self.nor_output.get())
        self.nor_output_edit.setToolTip('Suffix of the output FITS normalised spectra')        
        in_layout.addWidget(sfx_norm_label,3,0)
        in_layout.addWidget(self.nor_output_edit,3,1)

        inout_group.setLayout(in_layout)


        # Normalisation parameters
        norpar_group = QtWidgets.QGroupBox('Normalisation Parameters')

        # Normalisation parameters grid
        norpar_layout = QtWidgets.QGridLayout()

        # Spectra format (instrument and columns)
        instr_label = QtWidgets.QLabel("Instrument:")
        self.instruments_combo = QtWidgets.QComboBox()
        self.instruments_combo.setToolTip('Select an instrument to automatically define the columns/fields/hdus of the data.\nDEFAULT may be used for:\n- ASCII files with wavelength and flux (2 columns)\n- FITS monodimensional file with the wavelength from CRVAL1, CDELT1, NAXIS1\n- all HARPS/HARPS-N/SOPHIE s1d and e2ds data (old DRS)\n  and GIANO-B s1d data')        
        self.instruments_combo.addItems(self.option_values)
        self.instruments_combo.setCurrentText(self.option_instr.get())
        self.instruments_combo.currentTextChanged.connect(self.change_cols)

        wavecol_label = QtWidgets.QLabel("Wave:", alignment=Qt.AlignmentFlag.AlignRight)
        self.wavecol_spin = QtWidgets.QSpinBox(); self.wavecol_spin.setRange(0, 99); self.wavecol_spin.setValue(self.wavecol.get())
        self.wavecol_spin.setToolTip('Column/Field with wavelength (1 = 1st column)')
        fluxcol_label = QtWidgets.QLabel("Flux:", alignment=Qt.AlignmentFlag.AlignRight)
        self.fluxcol_spin = QtWidgets.QSpinBox(); self.fluxcol_spin.setRange(0, 99); self.fluxcol_spin.setValue(self.fluxcol.get())
        self.fluxcol_spin.setToolTip('Column/Field with flux (1 = 1st column)')
        snrcol_label = QtWidgets.QLabel("S/N:", alignment=Qt.AlignmentFlag.AlignRight)
        self.snrcol_spin = QtWidgets.QSpinBox(); self.snrcol_spin.setRange(0, 99); self.snrcol_spin.setValue(self.snrcol.get())
        self.snrcol_spin.setToolTip('Column/Field with S/N (1 = 1st column) - OPTIONAL')
        errcol_label = QtWidgets.QLabel("Errors:", alignment=Qt.AlignmentFlag.AlignRight)
        self.errcol_spin = QtWidgets.QSpinBox(); self.errcol_spin.setRange(0, 99); self.errcol_spin.setValue(self.errcol.get())
        self.errcol_spin.setToolTip('Column/Field with errors (1 = 1st column) - OPTIONAL')
        echcol_label = QtWidgets.QLabel("Orders:", alignment=Qt.AlignmentFlag.AlignRight)
        self.echcol_spin = QtWidgets.QSpinBox(); self.echcol_spin.setRange(0, 99); self.echcol_spin.setValue(self.echcol.get())
        self.echcol_spin.setToolTip('Column/Field with echelle orders (1 = 1st column) - OPTIONAL')

        
        norpar_layout.addWidget(instr_label,0,0)
        norpar_layout.addWidget(self.instruments_combo,0,1, 1,3)
        norpar_layout.addWidget(wavecol_label,1,0)
        norpar_layout.addWidget(self.wavecol_spin,1,1)
        norpar_layout.addWidget(fluxcol_label,1,2)
        norpar_layout.addWidget(self.fluxcol_spin,1,3)
        norpar_layout.addWidget(snrcol_label,1,4)
        norpar_layout.addWidget(self.snrcol_spin,1,5)
        norpar_layout.addWidget(errcol_label,1,6)
        norpar_layout.addWidget(self.errcol_spin,1,7)
        norpar_layout.addWidget(echcol_label,1,8)
        norpar_layout.addWidget(self.echcol_spin,1,9)        

        # Normalisation polynomials

        units_label = QtWidgets.QLabel("Wave units:")
        self.units_combo = QtWidgets.QComboBox(); self.units_combo.addItems(self.units_values); self.units_combo.setCurrentText(self.units_default.get()); self.units_combo.currentTextChanged.connect(self.change_units)
        self.units_combo.setToolTip('Select the wavelength unit')
        frame_label = QtWidgets.QLabel("Ref. frame:", alignment=Qt.AlignmentFlag.AlignRight)
        self.frame_combo = QtWidgets.QComboBox(); self.frame_combo.addItems(self.wave_values); self.frame_combo.setCurrentText(self.wave_frame.get()); self.frame_combo.currentTextChanged.connect(self.change_frame)
        self.frame_combo.setToolTip('Select the wavelength reference frame')
        deg_label = QtWidgets.QLabel("Poly degree:", alignment=Qt.AlignmentFlag.AlignRight)
        self.degree_spin = QtWidgets.QSpinBox(); self.degree_spin.setRange(0, 20); self.degree_spin.setValue(self.degree.get())
        self.degree_spin.setToolTip('Degree of polynomial for continuum fitting')
        ord_label = QtWidgets.QLabel("Subsets:", alignment=Qt.AlignmentFlag.AlignRight)
        self.n_ord_spin = QtWidgets.QSpinBox(); self.n_ord_spin.setRange(0, 99); self.n_ord_spin.setValue(self.n_ord.get())
        self.n_ord_spin.setToolTip('ONLY without echelle orders information:\ngive a number  of subsets to normalise independently')

        norpar_layout.addWidget(units_label,2,0)
        norpar_layout.addWidget(self.units_combo,2,1)
        norpar_layout.addWidget(frame_label,2,2)
        norpar_layout.addWidget(self.frame_combo,2,3)
        norpar_layout.addWidget(deg_label,2,4)
        norpar_layout.addWidget(self.degree_spin,2,5)
        norpar_layout.addWidget(ord_label,2,6)
        norpar_layout.addWidget(self.n_ord_spin,2,7)


 
        normalise_button = QtWidgets.QPushButton('Normalise spectra')
        normalise_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        normalise_button.clicked.connect(self.thread_normalise)
        
        norpar_layout.addWidget(normalise_button,3,0,1,10)
        
        norpar_group.setLayout(norpar_layout)
        
        #do_norm = QtWidgets.QHBoxLayout()
        #do_norm.addWidget(normalise_button)
        
        
        # Reset, abort, quit
        quit_group = QtWidgets.QGroupBox()

        # quit parameters grid
        quit_layout = QtWidgets.QGridLayout()
        
        reset_button = QtWidgets.QPushButton('Reset')
        reset_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        reset_button.clicked.connect(self.reset)
        reset_button.setToolTip('Reset all the fields to the default values\nClear the log and plot windows')
        
        abort_button = QtWidgets.QPushButton('ABORT')
        abort_button.setStyleSheet(f'QPushButton {{background-color: {self.color}; color: {self.back}; font: bold;}}')
        abort_button.clicked.connect(self.abort)
        abort_button.setToolTip('Abort the current process without exiting SHIVA')
        
        quit_button = QtWidgets.QPushButton('Quit')
        quit_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        quit_button.clicked.connect(self.quit_shiva)
        quit_button.setToolTip('Quit SHIVA, aborting any current process')
        
        quit_layout.addWidget(reset_button,0,0)
        quit_layout.addWidget(abort_button,0,1)
        quit_layout.addWidget(quit_button,0,2)
        
        quit_group.setLayout(quit_layout)
        

        #layout_all = QtWidgets.QVBoxLayout()
        layout.addWidget(inout_group)
        layout.addWidget(norpar_group)
        layout.addWidget(quit_group)
        
        #layout.addLayout(layout_all)
        #layout.addLayout(do_norm)
        layout.addStretch()
        return tab

    def create_profile_tab(self):
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        
        # Input/output parameters
        inout_group = QtWidgets.QGroupBox('Input/Output')
        
        # input grid
        in_layout = QtWidgets.QGridLayout()
        
        input_prf_label = QtWidgets.QLabel("Input folder:")        
        self.prf_indir_edit = QtWidgets.QLineEdit(self.prf_indir.get())
        self.prf_indir_edit.setToolTip('Directory with input normalised spectra')
        browse_dir = QtWidgets.QPushButton('Browse')
        browse_dir.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_dir.clicked.connect(self.prf_load_indir)
        
        in_layout.addWidget(input_prf_label,0,0)
        in_layout.addWidget(self.prf_indir_edit,0,1)
        in_layout.addWidget(browse_dir,0,2)

        # File/pattern folder
        spec_prf_label = QtWidgets.QLabel("File/Pattern:")        
        self.prf_spec_edit = QtWidgets.QLineEdit(self.prf_spec.get())
        self.prf_spec_edit.setToolTip('Select pattern OR single normalised spectrum')
        browse_file = QtWidgets.QPushButton('Browse')
        browse_file.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_file.clicked.connect(self.prf_load_file)
        in_layout.addWidget(spec_prf_label,1,0)
        in_layout.addWidget(self.prf_spec_edit,1,1)
        in_layout.addWidget(browse_file,1,2)

        # Output folder
        out_prf_label = QtWidgets.QLabel("Output folder:")
        self.prf_outdir_edit = QtWidgets.QLineEdit(self.nor_outdir.get())
        self.prf_outdir_edit.setToolTip('Select output folder')        
        browse_out = QtWidgets.QPushButton('Browse')
        browse_out.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_out.clicked.connect(self.prf_load_outdir)

        in_layout.addWidget(out_prf_label,2,0)
        in_layout.addWidget(self.prf_outdir_edit,2,1)
        in_layout.addWidget(browse_out,2,2)
        

        inout_group.setLayout(in_layout)
        
        # Single line extraction 
        ext_group = QtWidgets.QGroupBox('Single Line Extraction')
        # grid
        ext_layout = QtWidgets.QGridLayout()

        wave_label = QtWidgets.QLabel("Wavelength:", alignment=Qt.AlignmentFlag.AlignRight)
        self.ext_wave_spin = QtWidgets.QDoubleSpinBox(); self.ext_wave_spin.setRange(0.0, 1e7); self.ext_wave_spin.setValue(self.ext_wave.get())
        self.ext_wave_spin.setToolTip('Central wavelength to extract')
        rvmin_label = QtWidgets.QLabel("RV min:", alignment=Qt.AlignmentFlag.AlignRight)
        self.rvmin_ext_spin = QtWidgets.QDoubleSpinBox(); self.rvmin_ext_spin.setRange(-1e5, 1e5); self.rvmin_ext_spin.setValue(self.rvmin.get())
        self.rvmin_ext_spin.setToolTip('Minimum RV of the profile')
        rvmax_label = QtWidgets.QLabel("RV max:", alignment=Qt.AlignmentFlag.AlignRight)
        self.rvmax_ext_spin = QtWidgets.QDoubleSpinBox(); self.rvmax_ext_spin.setRange(-1e5, 1e5); self.rvmax_ext_spin.setValue(self.rvmax.get())
        self.rvmax_ext_spin.setToolTip('Maximum RV of the profile')
        rvstep_label = QtWidgets.QLabel("RV step:", alignment=Qt.AlignmentFlag.AlignRight)
        self.rvstep_ext_spin = QtWidgets.QDoubleSpinBox(); self.rvstep_ext_spin.setRange(0.001, 1e4); self.rvstep_ext_spin.setSingleStep(0.1); self.rvstep_ext_spin.setValue(self.rvstep.get())
        self.rvstep_ext_spin.setToolTip('RV step of the profile')
        
        self.rvmin_ext_spin.valueChanged.connect(self.change_rvmin)
        self.rvmax_ext_spin.valueChanged.connect(self.change_rvmax)
        self.rvstep_ext_spin.valueChanged.connect(self.change_rvstep)
        
        # Suffix
        sfx_label = QtWidgets.QLabel("Output suffix:")
        self.ext_output_edit = QtWidgets.QLineEdit(self.ext_output.get())
        self.ext_output_edit.setToolTip('Suffix of the output FITS profiles')

        extract_button = QtWidgets.QPushButton('Extract line')
        extract_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        extract_button.clicked.connect(self.thread_extract_line)

        ext_layout.addWidget(wave_label,0,0)
        ext_layout.addWidget(self.ext_wave_spin,0,1)
        ext_layout.addWidget(rvmin_label,0,2)
        ext_layout.addWidget(self.rvmin_ext_spin,0,3)
        ext_layout.addWidget(rvmax_label,0,4)
        ext_layout.addWidget(self.rvmax_ext_spin,0,5)
        ext_layout.addWidget(rvstep_label,0,6)
        ext_layout.addWidget(self.rvstep_ext_spin,0,7)
        ext_layout.addWidget(sfx_label,1,0)
        ext_layout.addWidget(self.ext_output_edit,1,1,1,7)
        
        ext_layout.addWidget(extract_button,2,0,1,8)
        ext_group.setLayout(ext_layout)
        

        # Mean line profile
        prf_group = QtWidgets.QGroupBox('Mean Line Profile Computation')
        
        # grid
        prf_layout = QtWidgets.QGridLayout()
        
        mask_label = QtWidgets.QLabel("Mask:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_edit = QtWidgets.QLineEdit(self.mask.get())
        self.mask_edit.setToolTip('Select mask:\nVALD file\n2-column ASCII file\nstandard FITS monodimensional spectrum)')
        browse_mask = QtWidgets.QPushButton('Browse')
        browse_mask.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_mask.clicked.connect(self.mask_load_file)        
        
        mask_invert_label = QtWidgets.QLabel("Invert:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_invert_box = QtWidgets.QCheckBox(); self.mask_invert_box.setChecked(bool(self.mask_invert.get()))
        self.mask_invert_box.setToolTip('Convert fluxes to depths (absorption mask/spectrum/model ONLY)')
        mask_spectrum_label = QtWidgets.QLabel("Spectrum:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_spectrum_box = QtWidgets.QCheckBox(); self.mask_spectrum_box.setChecked(bool(self.mask_spectrum.get()))
        self.mask_spectrum_box.setToolTip('Use a spectrum (with continuum) as mask -- CCF ONLY')

        mask_cwave_label = QtWidgets.QLabel("Wave col:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_wavecol_spin = QtWidgets.QSpinBox(); self.mask_wavecol_spin.setRange(0, 99); self.mask_wavecol_spin.setValue(self.mask_cwave.get())
        self.mask_wavecol_spin.setToolTip('Column/Field with wavelength (1 = 1st column)')
        mask_cflux_label = QtWidgets.QLabel("Flux col:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_fluxcol_spin = QtWidgets.QSpinBox(); self.mask_fluxcol_spin.setRange(0, 99); self.mask_fluxcol_spin.setValue(self.mask_cflux.get())
        self.mask_fluxcol_spin.setToolTip('Column/Field with wavelength (1 = 1st column)') 


        prf_layout.addWidget(mask_label,0,0) 
        prf_layout.addWidget(self.mask_edit,0,1,1,6) 
        prf_layout.addWidget(browse_mask,0,7)
                 
        prf_layout.addWidget(mask_invert_label,1,0) 
        prf_layout.addWidget(self.mask_invert_box,1,1) 
        prf_layout.addWidget(mask_spectrum_label,1,2)  
        prf_layout.addWidget(self.mask_spectrum_box,1,3)   
        prf_layout.addWidget(mask_cwave_label,1,4) 
        prf_layout.addWidget(self.mask_wavecol_spin,1,5) 
        prf_layout.addWidget(mask_cflux_label,1,6)  
        prf_layout.addWidget(self.mask_fluxcol_spin,1,7) 

        mask_units_label = QtWidgets.QLabel("Wave units:")
        self.mask_units_combo = QtWidgets.QComboBox(); self.mask_units_combo.addItems(self.mask_units_values); self.mask_units_combo.setCurrentText(self.mask_units_default.get()); self.mask_units_combo.currentTextChanged.connect(self.mask_change_units)
        self.mask_units_combo.setToolTip('Select the wavelength unit')
        mask_frame_label = QtWidgets.QLabel("Ref. frame:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_frame_combo = QtWidgets.QComboBox(); self.mask_frame_combo.addItems(self.mask_wave_values); self.mask_frame_combo.setCurrentText(self.mask_wave_frame.get()); self.mask_frame_combo.currentTextChanged.connect(self.mask_change_frame)
        self.mask_frame_combo.setToolTip('Select the wavelength reference frame')

        prf_layout.addWidget(mask_units_label,2,0) 
        prf_layout.addWidget(self.mask_units_combo,2,1,1,3) 
        prf_layout.addWidget(mask_frame_label,2,4) 
        prf_layout.addWidget(self.mask_frame_combo,2,5,1,3) 

        
        mask_dlow_label = QtWidgets.QLabel("Min. depth:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_dlow_spin = QtWidgets.QDoubleSpinBox(); self.mask_dlow_spin.setRange(-1000.0, 1.0); self.mask_dlow_spin.setSingleStep(0.05); self.mask_dlow_spin.setValue(self.mask_dlow.get())
        self.mask_dlow_spin.setToolTip('Minimum line depth of mask lines to be used')
        mask_dlup_label = QtWidgets.QLabel("Max. depth:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_dup_spin = QtWidgets.QDoubleSpinBox(); self.mask_dup_spin.setRange(0.01, 1000.0); self.mask_dup_spin.setSingleStep(0.05); self.mask_dup_spin.setValue(self.mask_dup.get())
        self.mask_dup_spin.setToolTip('Maximum line depth of mask lines to be used')
        mask_wmin_label = QtWidgets.QLabel("Min. wave:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_wmin_spin = QtWidgets.QDoubleSpinBox(); self.mask_wmin_spin.setRange(0.0, 1e6); self.mask_wmin_spin.setValue(self.mask_wmin.get())
        self.mask_wmin_spin.setToolTip('Minimum line wavelength of mask lines to be used')
        mask_wmax_label = QtWidgets.QLabel("Max. wave:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_wmax_spin = QtWidgets.QDoubleSpinBox(); self.mask_wmax_spin.setRange(0.0, 1e6); self.mask_wmax_spin.setValue(self.mask_wmax.get())
        self.mask_wmax_spin.setToolTip('Maximum line wavelength of mask lines to be used')


        prf_layout.addWidget(mask_dlow_label,3,0) 
        prf_layout.addWidget(self.mask_dlow_spin,3,1) 
        prf_layout.addWidget(mask_dlup_label,3,2) 
        prf_layout.addWidget(self.mask_dup_spin,3,3) 
        prf_layout.addWidget(mask_wmin_label,3,4)  
        prf_layout.addWidget(self.mask_wmin_spin,3,5) 
        prf_layout.addWidget(mask_wmax_label,3,6)  
        prf_layout.addWidget(self.mask_wmax_spin,3,7)

        
        mask_balmer_label = QtWidgets.QLabel("Exclude Balmer regions:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_balmer_box = QtWidgets.QCheckBox(); self.mask_balmer_box.setChecked(bool(self.mask_balmer.get()))
        self.mask_balmer_box.setToolTip('Balmer lines regions will not be used to compute the mean line profiles')
        mask_tell_label = QtWidgets.QLabel("Exclude telluric regions:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_tell_box = QtWidgets.QCheckBox(); self.mask_tell_box.setChecked(bool(self.mask_tell.get()))
        self.mask_tell_box.setToolTip('Regions heavily affected from telluric lines will not be used to compute the mean line profiles')
        prf_cosmic_label = QtWidgets.QLabel("Cosmics:", alignment=Qt.AlignmentFlag.AlignRight)
        self.cosmic_box = QtWidgets.QCheckBox(); self.cosmic_box.setChecked(bool(self.cosmic.get()))
        self.cosmic_box.setToolTip('Remove cosmics from spectra via sigma clipping prior to compute the mean line profiles')
        prf_clean_label = QtWidgets.QLabel("Smmoth spectra:", alignment=Qt.AlignmentFlag.AlignRight)
        self.clean_box = QtWidgets.QCheckBox(); self.clean_box.setChecked(bool(self.clean.get()))
        self.clean_box.setToolTip('Apply a smoothing spline to the spectra prior to compute the mean line profiles')
        prf_weight_label = QtWidgets.QLabel("S/N weigthed:", alignment=Qt.AlignmentFlag.AlignRight)
        self.ccfweight_box = QtWidgets.QCheckBox(); self.ccfweight_box.setChecked(bool(self.ccfweight.get()))
        self.ccfweight_box.setToolTip('CCF ONLY: use the normalised S/N values as weigths')

        prf_layout.addWidget(mask_balmer_label,4,0,1,3) 
        prf_layout.addWidget(self.mask_balmer_box,4,3) 
        prf_layout.addWidget(mask_tell_label,4,4,1,3) 
        prf_layout.addWidget(self.mask_tell_box,4,7) 
        
        prf_layout.addWidget(prf_cosmic_label,5,0)  
        prf_layout.addWidget(self.cosmic_box,5,1) 
        prf_layout.addWidget(prf_clean_label,5,2,1,2)  
        prf_layout.addWidget(self.clean_box,5,4)
        prf_layout.addWidget(prf_weight_label,5,5,1,2)  
        prf_layout.addWidget(self.ccfweight_box,5,7)

        
        mask_els_label = QtWidgets.QLabel("Select VALD elements:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_els_edit = QtWidgets.QLineEdit(self.mask_els.get())
        self.mask_els_edit.setToolTip('ONLY for VALD mask: use only selected elements, e.g. "Fe 1,Fe 2,H 1"')
        mask_noels_label = QtWidgets.QLabel("Exclude VALD elements:", alignment=Qt.AlignmentFlag.AlignRight)
        self.mask_noels_edit = QtWidgets.QLineEdit(self.mask_noels.get())
        self.mask_noels_edit.setToolTip('ONLY for VALD mask: exclude only selected elements, e.g. "Fe 1,Fe 2,H 1"')
        

        prf_layout.addWidget(mask_els_label,6,0,1,2) 
        prf_layout.addWidget(self.mask_els_edit,6,2,1,6) 
        prf_layout.addWidget(mask_noels_label,7,0,1,2) 
        prf_layout.addWidget(self.mask_noels_edit,7,2,1,6)         

        prf_rvmin_label = QtWidgets.QLabel("RV min:", alignment=Qt.AlignmentFlag.AlignRight)
        self.rvmin_prf_spin = QtWidgets.QDoubleSpinBox(); self.rvmin_prf_spin.setRange(-1e5, 1e5); self.rvmin_prf_spin.setValue(self.rvmin.get())
        self.rvmin_prf_spin.setToolTip('Minimum RV of the profile')
        prf_rvmax_label = QtWidgets.QLabel("RV max:", alignment=Qt.AlignmentFlag.AlignRight)
        self.rvmax_prf_spin = QtWidgets.QDoubleSpinBox(); self.rvmax_prf_spin.setRange(-1e5, 1e5); self.rvmax_prf_spin.setValue(self.rvmax.get())
        self.rvmax_prf_spin.setToolTip('Maximum RV of the profile')
        prf_rvstep_label = QtWidgets.QLabel("RV step:", alignment=Qt.AlignmentFlag.AlignRight)
        self.rvstep_prf_spin = QtWidgets.QDoubleSpinBox(); self.rvstep_prf_spin.setRange(0.001, 1e4); self.rvstep_prf_spin.setSingleStep(0.1); self.rvstep_prf_spin.setValue(self.rvstep.get())
        self.rvstep_prf_spin.setToolTip('RV step of the profile')
        
        self.rvmin_prf_spin.valueChanged.connect(self.change_rvmin)
        self.rvmax_prf_spin.valueChanged.connect(self.change_rvmax)
        self.rvstep_prf_spin.valueChanged.connect(self.change_rvstep)


        prf_layout.addWidget(prf_rvmin_label,8,0) 
        prf_layout.addWidget(self.rvmin_prf_spin,8,1) 
        prf_layout.addWidget(prf_rvmax_label,8,2) 
        prf_layout.addWidget(self.rvmax_prf_spin,8,3)  
        prf_layout.addWidget(prf_rvstep_label,8,4) 
        prf_layout.addWidget(self.rvstep_prf_spin,8,5)  
        
        
        # Suffix
        sfx_prf_label = QtWidgets.QLabel("Output suffix:")
        self.prf_output_edit = QtWidgets.QLineEdit(self.prf_output.get())
        self.prf_output_edit.setToolTip('Suffix of the output FITS mean line profiles')

        prf_layout.addWidget(sfx_prf_label,9,0) 
        prf_layout.addWidget(self.prf_output_edit,9,1,1,7)  

        # Which profile to run
                
        lsd_button = QtWidgets.QPushButton('Compute LSD')
        lsd_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        lsd_button.clicked.connect(self.thread_do_lsd_profile)

        ccf_button = QtWidgets.QPushButton('Compute CCF')
        ccf_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        ccf_button.clicked.connect(self.thread_do_ccf_profile)

        prf_layout.addWidget(lsd_button,10,0,1,4) 
        prf_layout.addWidget(ccf_button,10,4,1,4) 
        prf_group.setLayout(prf_layout)

        # Reset, abort, quit
        quit_group = QtWidgets.QGroupBox()

        # quit parameters grid
        quit_layout = QtWidgets.QGridLayout()
        
        reset_button = QtWidgets.QPushButton('Reset')
        reset_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        reset_button.clicked.connect(self.reset)
        reset_button.setToolTip('Reset all the fields to the default values\nClear the log and plot windows')
        
        abort_button = QtWidgets.QPushButton('ABORT')
        abort_button.setStyleSheet(f'QPushButton {{background-color: {self.color}; color: {self.back}; font: bold;}}')
        abort_button.clicked.connect(self.abort)
        abort_button.setToolTip('Abort the current process without exiting SHIVA')
        
        quit_button = QtWidgets.QPushButton('Quit')
        quit_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        quit_button.clicked.connect(self.quit_shiva)
        quit_button.setToolTip('Quit SHIVA, aborting any current process')
        
        quit_layout.addWidget(reset_button,0,0)
        quit_layout.addWidget(abort_button,0,1)
        quit_layout.addWidget(quit_button,0,2)
        
        quit_group.setLayout(quit_layout)


        #layout_all = QtWidgets.QVBoxLayout()
        layout.addWidget(inout_group)
        layout.addWidget(ext_group)
        layout.addWidget(prf_group)
        layout.addWidget(quit_group)
        
        #layout.addLayout(layout_all)

        layout.addStretch()
        
        return tab

    def create_analysis_tab(self):
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        
        # Read profiles
        file_group = QtWidgets.QGroupBox('Input/Output')
                
        file_layout = QtWidgets.QGridLayout()

        input_la_label = QtWidgets.QLabel("Input folder:")        
        self.la_indir_edit = QtWidgets.QLineEdit(self.la_indir.get())
        self.la_indir_edit.setToolTip('Directory with input line profiles')
        browse_dir = QtWidgets.QPushButton('Browse')
        browse_dir.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_dir.clicked.connect(self.la_load_indir)
        
        file_layout.addWidget(input_la_label,0,0)
        file_layout.addWidget(self.la_indir_edit,0,1)
        file_layout.addWidget(browse_dir,0,2)

        # File/pattern folder
        spec_la_label = QtWidgets.QLabel("File/Pattern:")        
        self.la_spec_edit = QtWidgets.QLineEdit(self.la_spec.get())
        self.la_spec_edit.setToolTip('Select pattern OR single line profile')
        browse_file = QtWidgets.QPushButton('Browse')
        browse_file.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_file.clicked.connect(self.la_load_file)
        
        file_layout.addWidget(spec_la_label,1,0)
        file_layout.addWidget(self.la_spec_edit,1,1)
        file_layout.addWidget(browse_file,1,2)

        # Output folder
        out_la_label = QtWidgets.QLabel("Output folder:")
        self.la_outdir_edit = QtWidgets.QLineEdit(self.la_outdir.get())
        self.la_outdir_edit.setToolTip('Select output folder')        
        browse_out = QtWidgets.QPushButton('Browse')
        browse_out.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_out.clicked.connect(self.la_load_outdir)

        file_layout.addWidget(out_la_label,2,0)
        file_layout.addWidget(self.la_outdir_edit,2,1)
        file_layout.addWidget(browse_out,2,2)
        
        file_group.setLayout(file_layout)

        # Normalise profiles        
        norm_group = QtWidgets.QGroupBox('Profile Normalisation')
        
        norm_layout = QtWidgets.QGridLayout()

        la_low_label = QtWidgets.QLabel("RV lower line limit:")                
        self.limitlow_spin = QtWidgets.QDoubleSpinBox(); self.limitlow_spin.setRange(-1e5, 1e5); self.limitlow_spin.setValue(self.limitlow.get())
        self.limitlow_spin.setToolTip('Define the lower limit of the line, for continuum normalisation')
        la_up_label = QtWidgets.QLabel("RV upper line limit:")                
        self.limitup_spin = QtWidgets.QDoubleSpinBox(); self.limitup_spin.setRange(-1e5, 1e5); self.limitup_spin.setValue(self.limitup.get())
        self.limitup_spin.setToolTip('Define the upper limit of the line, for continuum normalisation')
        
        self.limitlow_spin.valueChanged.connect(self.la_change_rvmin)
        self.limitup_spin.valueChanged.connect(self.la_change_rvmax)
        
        la_std_label = QtWidgets.QLabel("Time series StDev:")                
        self.std_box = QtWidgets.QCheckBox(); self.std_box.setChecked(bool(self.std.get()))
        self.std_box.setToolTip('Compute the standard deviation of the line profiles from their average, and save the data')

        la_norprf_sfx_label = QtWidgets.QLabel("Output suffix:")                
        self.norprf_output_edit = QtWidgets.QLineEdit(self.norprf_output.get())
        self.norprf_output_edit.setToolTip('Suffix of the output FITS mean line profiles')
        
        norm_button = QtWidgets.QPushButton('Normalise profiles')
        norm_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        norm_button.clicked.connect(self.thread_norm_profile)
        
        norm_layout.addWidget(la_low_label,0,0,1,2)
        norm_layout.addWidget(self.limitlow_spin,0,2)
        norm_layout.addWidget(la_up_label,0,3,1,2)
        norm_layout.addWidget(self.limitup_spin,0,5)
        norm_layout.addWidget(la_std_label,0,6)
        norm_layout.addWidget(self.std_box,0,7)

        norm_layout.addWidget(la_norprf_sfx_label,1,0,1,2) 
        norm_layout.addWidget(self.norprf_output_edit,1,2,1,6)  
        
        norm_layout.addWidget(norm_button,2,0,1,8)
        norm_group.setLayout(norm_layout)

        # Fitting line profiles
        fit_group = QtWidgets.QGroupBox('Profile Fitting')
        
        fit_layout = QtWidgets.QGridLayout()
        #fit_layout = QtWidgets.QFormLayout()
        
        # Separate components
        la_first_label = QtWidgets.QLabel("First:")
        self.la_first_combo = QtWidgets.QComboBox(); self.la_first_combo.addItems(self.option_la_fit); self.la_first_combo.setCurrentText(self.option_fit.get()); self.la_first_combo.currentTextChanged.connect(self.la_change_fit)
        self.la_first_combo.setToolTip('Select the fitting function for the first (or only) component')

        la_firstRV_label = QtWidgets.QLabel("RV guess:")
        self.rv0_spin = QtWidgets.QDoubleSpinBox(); self.rv0_spin.setRange(-1e5, 1e5); self.rv0_spin.setValue(self.rv0.get())
        self.rv0_spin.setToolTip('RV guess in km/s for the first (or only) component')
        
        la_firstwidth_label = QtWidgets.QLabel("Width guess:")
        self.width_spin = QtWidgets.QDoubleSpinBox(); self.width_spin.setRange(0, 1e5); self.width_spin.setValue(self.width.get())
        self.width_spin.setToolTip('Width guess in km/s for the first (or only) component')
        
        la_firstLD_label = QtWidgets.QLabel("Linear LD:")        
        self.ld_spin = QtWidgets.QDoubleSpinBox(); self.ld_spin.setRange(0, 1.0); self.ld_spin.setSingleStep(0.01); self.ld_spin.setValue(self.ld.get())
        self.ld_spin.setToolTip('Linear limb darkening (fixed value) for the first (or only) component')
        
        fit_layout.addWidget(la_first_label,0,0)
        fit_layout.addWidget(self.la_first_combo,0,1)
        fit_layout.addWidget(la_firstRV_label,0,2)
        fit_layout.addWidget(self.rv0_spin,0,3)
        fit_layout.addWidget(la_firstwidth_label,0,4)
        fit_layout.addWidget(self.width_spin,0,5)
        fit_layout.addWidget(la_firstLD_label,0,6)
        fit_layout.addWidget(self.ld_spin,0,7)
        

        la_second_label = QtWidgets.QLabel("Second:")
        self.la_second_combo = QtWidgets.QComboBox(); self.la_second_combo.addItems(self.option_la_fits); self.la_second_combo.setCurrentText(self.option_fit2.get()); self.la_second_combo.currentTextChanged.connect(self.la_change_fit2)
        self.la_second_combo.setToolTip('Select the fitting function for the second component (OPTIONAL)')
        
        la_secondRV_label = QtWidgets.QLabel("RV guess:")
        self.rv02_spin = QtWidgets.QDoubleSpinBox(); self.rv02_spin.setRange(-1e5, 1e5); self.rv02_spin.setValue(self.rv02.get())
        self.rv02_spin.setToolTip('RV guess in km/s for the second component (if any)')
        
        la_secondwidth_label = QtWidgets.QLabel("Width guess:")
        self.width2_spin = QtWidgets.QDoubleSpinBox(); self.width2_spin.setRange(0, 1e5); self.width2_spin.setValue(self.width2.get())
        self.width2_spin.setToolTip('Width guess in km/s for the second component (if any)')
        
        la_secondLD_label = QtWidgets.QLabel("Linear LD:")        
        self.ld2_spin = QtWidgets.QDoubleSpinBox(); self.ld2_spin.setRange(0, 1.0); self.ld2_spin.setSingleStep(0.01); self.ld2_spin.setValue(self.ld2.get())
        self.ld2_spin.setToolTip('Linear limb darkening (fixed value) for the second component (if any)')
        
        fit_layout.addWidget(la_second_label,1,0)
        fit_layout.addWidget(self.la_second_combo,1,1)
        fit_layout.addWidget(la_secondRV_label,1,2)
        fit_layout.addWidget(self.rv02_spin,1,3)
        fit_layout.addWidget(la_secondwidth_label,1,4)
        fit_layout.addWidget(self.width2_spin,1,5)
        fit_layout.addWidget(la_secondLD_label,1,6)
        fit_layout.addWidget(self.ld2_spin,1,7)

        
        la_third_label = QtWidgets.QLabel("Third:")
        self.la_third_combo = QtWidgets.QComboBox(); self.la_third_combo.addItems(self.option_la_fits); self.la_third_combo.setCurrentText(self.option_fit3.get()); self.la_third_combo.currentTextChanged.connect(self.la_change_fit3)
        self.la_third_combo.setToolTip('Select the fitting function for the third component (OPTIONAL)')

        la_thirdRV_label = QtWidgets.QLabel("RV guess:")
        self.rv03_spin = QtWidgets.QDoubleSpinBox(); self.rv03_spin.setRange(-1e5, 1e5); self.rv03_spin.setValue(self.rv03.get())
        self.rv03_spin.setToolTip('RV guess in km/s for the third component (if any)')
        
        la_thirdwidth_label = QtWidgets.QLabel("Width guess:")
        self.width3_spin = QtWidgets.QDoubleSpinBox(); self.width3_spin.setRange(0, 1e5); self.width3_spin.setValue(self.width3.get())
        self.width3_spin.setToolTip('Width guess in km/s for the third component (if any)')
        
        la_thirdLD_label = QtWidgets.QLabel("Linear LD:")        
        self.ld3_spin = QtWidgets.QDoubleSpinBox(); self.ld3_spin.setRange(0, 1.0); self.ld3_spin.setSingleStep(0.01); self.ld3_spin.setValue(self.ld3.get())
        self.ld3_spin.setToolTip('Linear limb darkening (fixed value) for the third component (if any)')
        
        fit_layout.addWidget(la_third_label,2,0)
        fit_layout.addWidget(self.la_third_combo,2,1)
        fit_layout.addWidget(la_thirdRV_label,2,2)
        fit_layout.addWidget(self.rv03_spin,2,3)
        fit_layout.addWidget(la_thirdwidth_label,2,4)
        fit_layout.addWidget(self.width3_spin,2,5)
        fit_layout.addWidget(la_thirdLD_label,2,6)
        fit_layout.addWidget(self.ld3_spin,2,7)


        la_errors_label = QtWidgets.QLabel("Use errors:")
        self.fit_errs_box = QtWidgets.QCheckBox(); self.fit_errs_box.setChecked(bool(self.fit_errs.get()))
        self.fit_errs_box.setToolTip('Use the errors of the data when fitting')
        
        la_res_label = QtWidgets.QLabel("Resolution:")
        self.res_spin = QtWidgets.QSpinBox(); self.res_spin.setRange(0, 1000000); self.res_spin.setSingleStep(1000); self.res_spin.setValue(self.resolution.get())
        self.res_spin.setToolTip('Instrumental resolution: if given, it will be considered for all fitting functions excluding rotational function and Supergaussian')

        fit_layout.addWidget(la_errors_label,3,0)
        fit_layout.addWidget(self.fit_errs_box,3,1)
        fit_layout.addWidget(la_res_label,3,2)
        fit_layout.addWidget(self.res_spin,3,3)


        la_fit_sfx_label = QtWidgets.QLabel("Save fit report:", alignment=Qt.AlignmentFlag.AlignRight)
        self.save_report_box = QtWidgets.QCheckBox(); self.save_report_box.setChecked(bool(self.save_report.get()))
        self.save_report_box.setToolTip('Save fit report as a text file')
        
        fit_button = QtWidgets.QPushButton('Fit line profiles')
        fit_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        fit_button.clicked.connect(self.thread_fit_profile)
        
        fit_layout.addWidget(la_fit_sfx_label,3,4,1,3)
        fit_layout.addWidget(self.save_report_box,3,7)
        fit_layout.addWidget(fit_button,4,0,1,8)        

        fit_group.setLayout(fit_layout)


        # Derive the line moments/bisector/FT from the profiles
        other_group = QtWidgets.QGroupBox('Line analysis: Moments, Bisector, Fourier Transform')
        other_layout = QtWidgets.QGridLayout()

        la_limits_label = QtWidgets.QLabel("Limits:")
        self.la_limits_combo = QtWidgets.QComboBox(); self.la_limits_combo.addItems(self.option_la_values); self.la_limits_combo.setCurrentText(self.option_limits.get()); self.la_limits_combo.currentTextChanged.connect(self.option_la_change)
        self.la_limits_combo.setToolTip('Define the line limits by the prior fit or manually')
        
        la_low_lim_label = QtWidgets.QLabel("RV lower line limit:")
        self.la_limitlow_spin = QtWidgets.QDoubleSpinBox(); self.la_limitlow_spin.setRange(-1e5, 1e5); self.la_limitlow_spin.setValue(self.limitlow.get())
        self.la_limitlow_spin.setToolTip('Manually define the lower limit of the line, for line analysis - ONLY if Manual limits are selected')
        la_up_lim_label = QtWidgets.QLabel("RV upper line limit:")                
        self.la_limitup_spin = QtWidgets.QDoubleSpinBox(); self.la_limitup_spin.setRange(-1e5, 1e5); self.la_limitup_spin.setValue(self.limitup.get())
        self.la_limitup_spin.setToolTip('Manually define the upper limit of the line, for line analysis - ONLY if Manual limits are selected')

        self.la_limitlow_spin.valueChanged.connect(self.la_change_rvmin)
        self.la_limitup_spin.valueChanged.connect(self.la_change_rvmax)

        other_layout.addWidget(la_limits_label,0,0)
        other_layout.addWidget(self.la_limits_combo,0,1)
        other_layout.addWidget(la_low_lim_label,0,2)
        other_layout.addWidget(self.la_limitlow_spin,0,3)
        other_layout.addWidget(la_up_lim_label,0,4)
        other_layout.addWidget(self.la_limitup_spin,0,5)

        
        la_fou_sfx_label = QtWidgets.QLabel("Fourier suffix:")                
        self.fou_output_edit = QtWidgets.QLineEdit(self.fou_output.get())
        self.fou_output_edit.setToolTip('Suffix of the output FITS Fourier results')

        mom_button = QtWidgets.QPushButton('Moments')
        mom_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        mom_button.clicked.connect(self.thread_moments)
        
        bis_button = QtWidgets.QPushButton('Bisector')
        bis_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        bis_button.clicked.connect(self.thread_bisector)
        
        fou_button = QtWidgets.QPushButton('Fourier')
        fou_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        fou_button.clicked.connect(self.thread_fourier)

        other_layout.addWidget(la_fou_sfx_label,1,0)
        other_layout.addWidget(self.fou_output_edit,1,1,1,5)        

        other_layout.addWidget(mom_button,2,0,1,2)
        
        other_layout.addWidget(bis_button,2,2,1,2)
        
        other_layout.addWidget(fou_button,2,4,1,2)

        other_group.setLayout(other_layout)

        # Reset, abort, quit
        quit_group = QtWidgets.QGroupBox()

        # quit parameters grid
        quit_layout = QtWidgets.QGridLayout()
        
        reset_button = QtWidgets.QPushButton('Reset')
        reset_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        reset_button.clicked.connect(self.reset)
        reset_button.setToolTip('Reset all the fields to the default values\nClear the log and plot windows')
        
        abort_button = QtWidgets.QPushButton('ABORT')
        abort_button.setStyleSheet(f'QPushButton {{background-color: {self.color}; color: {self.back}; font: bold;}}')
        abort_button.clicked.connect(self.abort)
        abort_button.setToolTip('Abort the current process without exiting SHIVA')
        
        quit_button = QtWidgets.QPushButton('Quit')
        quit_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        quit_button.clicked.connect(self.quit_shiva)
        quit_button.setToolTip('Quit SHIVA, aborting any current process')
        
        quit_layout.addWidget(reset_button,0,0)
        quit_layout.addWidget(abort_button,0,1)
        quit_layout.addWidget(quit_button,0,2)
        
        quit_group.setLayout(quit_layout)


        layout.addWidget(file_group)
        layout.addWidget(norm_group)
        layout.addWidget(fit_group)
        layout.addWidget(other_group)
        layout.addWidget(quit_group)
        
        layout.addStretch()
        return tab


    def create_timeseries_tab(self):
        tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(tab)
        
        # Read profiles
        file_group = QtWidgets.QGroupBox('Input/Output')
                
        file_layout = QtWidgets.QGridLayout()

        input_ts_label = QtWidgets.QLabel("Input folder:")        
        self.ts_indir_edit = QtWidgets.QLineEdit(self.ts_indir.get())
        self.ts_indir_edit.setToolTip('Directory with input analyzed line profiles')
        browse_dir = QtWidgets.QPushButton('Browse')
        browse_dir.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_dir.clicked.connect(self.ts_load_indir)
        
        file_layout.addWidget(input_ts_label,0,0)
        file_layout.addWidget(self.ts_indir_edit,0,1)
        file_layout.addWidget(browse_dir,0,2)

        # File/pattern folder
        spec_ts_label = QtWidgets.QLabel("File/Pattern:")        
        self.ts_spec_edit = QtWidgets.QLineEdit(self.ts_spec.get())
        self.ts_spec_edit.setToolTip('Select pattern OR single line profile')
        browse_file = QtWidgets.QPushButton('Browse')
        browse_file.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_file.clicked.connect(self.ts_load_file)
        
        file_layout.addWidget(spec_ts_label,1,0)
        file_layout.addWidget(self.ts_spec_edit,1,1)
        file_layout.addWidget(browse_file,1,2)

        # Output folder
        out_ts_label = QtWidgets.QLabel("Output folder:")
        self.ts_outdir_edit = QtWidgets.QLineEdit(self.ts_outdir.get())
        self.ts_outdir_edit.setToolTip('Select output folder')        
        browse_out = QtWidgets.QPushButton('Browse')
        browse_out.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        browse_out.clicked.connect(self.ts_load_outdir)

        file_layout.addWidget(out_ts_label,2,0)
        file_layout.addWidget(self.ts_outdir_edit,2,1)
        file_layout.addWidget(browse_out,2,2)
        
        file_group.setLayout(file_layout)

        # Plot time series        
        plot_group = QtWidgets.QGroupBox('Plot time series')
        
        plot_layout = QtWidgets.QGridLayout()

        ts_xdata_label = QtWidgets.QLabel("X data:")    
        
        self.ts_xdata_combo = QtWidgets.QComboBox(); self.ts_xdata_combo.addItems(self.option_ts_plot); self.ts_xdata_combo.setCurrentText(self.ts_xplot.get()); self.ts_xdata_combo.currentTextChanged.connect(self.ts_xdata_change)
        self.ts_xdata_combo.setToolTip('Data to be plotted on the x-axis')
        
        ts_ydata_label = QtWidgets.QLabel("Y data:")    
        
        self.ts_ydata_combo = QtWidgets.QComboBox(); self.ts_ydata_combo.addItems(self.option_ts_plot); self.ts_ydata_combo.setCurrentText(self.ts_yplot.get()); self.ts_ydata_combo.currentTextChanged.connect(self.ts_xdata_change)
        self.ts_ydata_combo.setToolTip('Data to be plotted on the y-axis')
                
        ts_error_label = QtWidgets.QLabel("Plot errors:")
        self.ts_error_box = QtWidgets.QCheckBox(); self.ts_error_box.setChecked(bool(self.ts_error.get()))
        self.ts_error_box.setToolTip('Plot with errorbars')
        
        ts_comp_label = QtWidgets.QLabel("Fit component:")
        self.ts_comp_combo = QtWidgets.QComboBox(); self.ts_comp_combo.addItems(self.option_ts_comp); self.ts_comp_combo.setCurrentText(self.ts_comp.get()); self.ts_comp_combo.currentTextChanged.connect(self.ts_comp_change)
        self.ts_comp_combo.setToolTip('Fitting component - used ONLY if the data are from the fitting function: RV, EW, width')
        
        ts_save_label = QtWidgets.QLabel("Save time series:")
        self.ts_save_box = QtWidgets.QCheckBox(); self.ts_save_box.setChecked(bool(self.ts_save.get()))
        self.ts_save_box.setToolTip('Save the data as a txt file with 4 columns:\nx_data, x_err, y_data, y_err')
                
        plot_button = QtWidgets.QPushButton('Plot time series')
        plot_button.setStyleSheet(f'QPushButton {{background-color: {self.runbutton};}}')
        plot_button.clicked.connect(self.thread_plot_ts)
        
        plot_layout.addWidget(ts_xdata_label,0,0)
        plot_layout.addWidget(self.ts_xdata_combo,0,1)
        plot_layout.addWidget(ts_ydata_label,0,2)
        plot_layout.addWidget(self.ts_ydata_combo,0,3)
        plot_layout.addWidget(ts_error_label,0,4)
        plot_layout.addWidget(self.ts_error_box,0,5)
        plot_layout.addWidget(ts_comp_label,0,6)
        plot_layout.addWidget(self.ts_comp_combo,0,7)
        plot_layout.addWidget(ts_save_label,0,8) 
        plot_layout.addWidget(self.ts_save_box,0,9)  
        
        plot_layout.addWidget(plot_button,1,0,1,10)
        plot_group.setLayout(plot_layout)

        # Reset, abort, quit
        quit_group = QtWidgets.QGroupBox()

        # quit parameters grid
        quit_layout = QtWidgets.QGridLayout()
        
        reset_button = QtWidgets.QPushButton('Reset')
        reset_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        reset_button.clicked.connect(self.reset)
        reset_button.setToolTip('Reset all the fields to the default values\nClear the log and plot windows')
        
        abort_button = QtWidgets.QPushButton('ABORT')
        abort_button.setStyleSheet(f'QPushButton {{background-color: {self.color}; color: {self.back}; font: bold;}}')
        abort_button.clicked.connect(self.abort)
        abort_button.setToolTip('Abort the current process without exiting SHIVA')
        
        quit_button = QtWidgets.QPushButton('Quit')
        quit_button.setStyleSheet(f'QPushButton {{background-color: {self.back}; color: {self.color};}}')
        quit_button.clicked.connect(self.quit_shiva)
        quit_button.setToolTip('Quit SHIVA, aborting any current process')
        
        quit_layout.addWidget(reset_button,0,0)
        quit_layout.addWidget(abort_button,0,1)
        quit_layout.addWidget(quit_button,0,2)
        
        quit_group.setLayout(quit_layout)


        layout.addWidget(file_group)
        layout.addWidget(plot_group)
        layout.addWidget(quit_group)
        
        layout.addStretch()
        return tab


    def create_plot_group(self):
        group = QtWidgets.QGroupBox('Plots')
        layout = QtWidgets.QVBoxLayout(group)
        self.fig, self.ax = plt.subplots(figsize=(5, 3), dpi=100)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas)
        layout.addWidget(self.canvas)
        layout.addWidget(self.toolbar)
        self.plot_logo()
        return group

    def create_log_group(self):
        group = QtWidgets.QGroupBox('Messages')
        #log_layout = QtWidgets.QGridLayout()
        
        layout = QtWidgets.QVBoxLayout(group)
        self.log_text = QtWidgets.QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setPlainText(self.messages.get())
        
        log_reset_button = QtWidgets.QPushButton('Reset Log')
        log_reset_button.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        log_reset_button.clicked.connect(self.reset_log)
        log_reset_button.setToolTip('Clear the log window')
        
        save_button = QtWidgets.QPushButton('Save Log')
        save_button.setStyleSheet(f'QPushButton {{background-color: {self.browse};}}')    
        save_button.clicked.connect(self.save_log)
        save_button.setToolTip('Save the current log as a text file')
        
        layout.addWidget(self.log_text)
        layout.addWidget(log_reset_button)
        layout.addWidget(save_button)
        return group

    def reset(self):
        # reset all fields, clear log and plot
        self.set_entries()
        self.log_text.clear()
        self.log_text.setPlainText(self.messages.get())
        self.plot_logo()
        return


    def abort(self):
        # end all processes
        if self.thread_running.get():
            self.abort_value.set(True)
            self.update_log("\n########################\nAborting the processes\nWaiting for the right moment...\n########################\n", timestamp=False)
        else:
            self.update_log("\n########################\nNo process is running.\n########################\n", timestamp=False)        
        return

        
    def quit_shiva(self):
        # quit the GUI
        self.abort()
        sys.exit()
        return


    def append_log(self, message):
        self.log_text.append(message)
        self.log_text.verticalScrollBar().setValue(self.log_text.verticalScrollBar().maximum())

    def update_log(self, message, timestamp=True):
        if timestamp:
            now = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
            formatted = ':\n'.join((now, message))
        else:
            formatted = message
        self.log_signal.emit(formatted)

    def plot_logo(self):
        self.ax.clear()
        try:
            img = plt.imread(self.sp_logo)
            self.ax.imshow(img)
            self.ax.axis('off')
        except FileNotFoundError:
            self.ax.plot([], [])
        self.canvas.draw()

    def update_figure(self, xvalues, yvalues, x_err=None, y_err=None, x_adds=None, y_adds=None, y_res=None, ymin=None, ymax=None, limits=(None, None), logscale=False):
        self.ax.clear()
        try:
            self.ax1.clear()
        except AttributeError:
            pass
        if y_res is None:
            if np.logical_and(x_err is None, y_err is None):
                self.ax.plot(xvalues, yvalues)
            else:
                self.ax.errorbar(xvalues, yvalues, xerr=x_err, yerr=y_err, fmt='o')
        else:
            if np.logical_and(x_err is None, y_err is None):
                self.ax.plot(xvalues, yvalues)
            else:
                self.ax.errorbar(xvalues, yvalues, xerr=x_err, yerr=y_err, fmt='o')
            self.ax.plot(xvalues, y_res)
        if y_adds is not None:
            if x_adds is not None:
                xvalues = x_adds
            for y_add in y_adds:
                self.ax.plot(xvalues, y_add)
        if limits is not None:
            for limit in limits:
                if limit is not None:
                    self.ax.axvline(limit, color='gray', linestyle='--')
        if ymin is not None or ymax is not None:
            self.ax.set_ylim(bottom=ymin, top=ymax)
        if logscale:
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')
        #self.ax.legend(loc='best')
        self.canvas.draw()

    def find_file(self, filepattern, fileplace=None):
        if os.path.isfile(filepattern):
            return [filepattern]
        wildcards = ''.join(('*', filepattern))
        if fileplace is None:
            fileplace = self.basedir
        specs = glob.glob(os.path.join(fileplace, filepattern))
        return sorted(specs)

    def save_fits(self, data, header, oldname, outdir, suffix=None):
        columns = []
        for key in data:
            if isinstance(data[key], np.ndarray):
                columns.append(fits.Column(name=key.upper(), format='D', array=data[key]))
            else:
                columns.append(fits.Column(name=key.upper(), format='D', array=np.ones(2) * data[key]))
        tbhdu = fits.BinTableHDU.from_columns(columns)
        prihdu = fits.PrimaryHDU(data=None, header=header)
        hdulist = fits.HDUList([prihdu, tbhdu])
        if not os.path.isdir(outdir):
            os.makedirs(outdir, exist_ok=True)
        if suffix:
            basename = os.path.splitext(os.path.basename(oldname))[0]
            newname = os.path.join(outdir, f"{basename}{suffix}")
        else:
            newname = oldname
        hdulist.writeto(newname, overwrite=True, output_verify='ignore')

    def insert_key(self, hea, keyword, value, comment):
        value = np.nan_to_num(value, nan=0.0, posinf=0.0, neginf=0.0)
        hea[keyword] = (value, comment)
        return hea

    def update_fits(self, fitsname, hea):
        with fits.open(fitsname, mode='update') as hdu:
            for entry in hea:
                try:
                    hdu[0].header[entry] = hea[entry]
                except ValueError:
                    pass
            hdu.flush()

    def sync_ui_values(self):
        self.nor_indir.set(self.nor_indir_edit.text())
        self.nor_spec.set(self.nor_spec_edit.text())
        self.option_instr.set(self.instruments_combo.currentText())
        self.wavecol.set(self.wavecol_spin.value())
        self.fluxcol.set(self.fluxcol_spin.value())
        self.snrcol.set(self.snrcol_spin.value())
        self.errcol.set(self.errcol_spin.value())
        self.echcol.set(self.echcol_spin.value())
        self.units_default.set(self.units_combo.currentText())
        self.wave_frame.set(self.frame_combo.currentText())
        self.degree.set(self.degree_spin.value())
        self.n_ord.set(self.n_ord_spin.value())
        self.nor_outdir.set(self.nor_outdir_edit.text())
        self.nor_output.set(self.nor_output_edit.text())
        self.prf_indir.set(self.prf_indir_edit.text())
        self.prf_spec.set(self.prf_spec_edit.text())
        self.prf_outdir.set(self.prf_outdir_edit.text())
        self.mask.set(self.mask_edit.text())
        self.mask_invert.set(int(self.mask_invert_box.isChecked()))
        self.mask_spectrum.set(int(self.mask_spectrum_box.isChecked()))
        self.mask_cwave.set(self.mask_wavecol_spin.value())
        self.mask_cflux.set(self.mask_fluxcol_spin.value())
        self.mask_units_default.set(self.mask_units_combo.currentText())
        self.mask_wave_frame.set(self.mask_frame_combo.currentText())
        self.mask_dlow.set(self.mask_dlow_spin.value())
        self.mask_dup.set(self.mask_dup_spin.value())
        self.mask_wmin.set(self.mask_wmin_spin.value())
        self.mask_wmax.set(self.mask_wmax_spin.value())
        self.mask_els.set(self.mask_els_edit.text())
        self.mask_noels.set(self.mask_noels_edit.text())
        self.mask_balmer.set(int(self.mask_balmer_box.isChecked()))
        self.mask_tell.set(int(self.mask_tell_box.isChecked()))
        self.cosmic.set(int(self.cosmic_box.isChecked()))
        self.clean.set(int(self.clean_box.isChecked()))
        self.ccfweight.set(int(self.ccfweight_box.isChecked()))
        self.ext_wave.set(self.ext_wave_spin.value())
        self.rvmin.set(self.rvmin_ext_spin.value())
        self.rvmax.set(self.rvmax_ext_spin.value())
        self.rvstep.set(self.rvstep_ext_spin.value())
        self.ext_output.set(self.ext_output_edit.text())
        self.prf_output.set(self.prf_output_edit.text())
        self.la_indir.set(self.la_indir_edit.text())
        self.la_spec.set(self.la_spec_edit.text())
        self.la_outdir.set(self.la_outdir_edit.text())
        self.limitlow.set(self.limitlow_spin.value())
        self.limitup.set(self.limitup_spin.value())
        self.std.set(int(self.std_box.isChecked()))
        self.norprf_output.set(self.norprf_output_edit.text())
        self.chosen_fit.set(self.la_first_combo.currentText())
        self.chosen_fit2.set(self.la_second_combo.currentText())
        self.chosen_fit3.set(self.la_third_combo.currentText())
        self.rv0.set(self.rv0_spin.value())
        self.width.set(self.width_spin.value())
        self.ld.set(self.ld_spin.value())
        self.fit_errs.set(int(self.fit_errs_box.isChecked()))
        self.save_report.set(int(self.save_report_box.isChecked()))
        #self.la_limitlow.set(self.la_limitlow_spin.value())
        #self.la_limitup.set(self.la_limitup_spin.value())
        self.fou_output.set(self.fou_output_edit.text())
        self.ts_indir.set(self.ts_indir_edit.text())
        self.ts_spec.set(self.ts_spec_edit.text())
        self.ts_outdir.set(self.ts_outdir_edit.text())
        self.ts_xplot.set(self.ts_xdata_combo.currentText())
        self.ts_yplot.set(self.ts_ydata_combo.currentText())
        self.ts_comp.set(self.ts_comp_combo.currentText())        
        self.ts_error.set(int(self.ts_error_box.isChecked()))
        self.ts_save.set(int(self.ts_save_box.isChecked()))

    def change_cols(self, selection):
        self.option_instr.set(selection)
        if selection == 'ESPRESSO S1D':
            self.wavecol.set(1); self.fluxcol.set(3); self.snrcol.set(0); self.errcol.set(4); self.echcol.set(0); self.spec_unit.set('a')
        elif selection == 'ESPRESSO S2D':
            self.wavecol.set(4); self.fluxcol.set(1); self.snrcol.set(0); self.errcol.set(2); self.echcol.set(0); self.spec_unit.set('a')
        elif selection == 'GIANO-B MS1D':
            self.wavecol.set(2); self.fluxcol.set(3); self.snrcol.set(4); self.errcol.set(0); self.echcol.set(1); self.spec_unit.set('n')
        elif selection == 'CARMENES':
            self.wavecol.set(4); self.fluxcol.set(1); self.snrcol.set(0); self.errcol.set(3); self.echcol.set(0); self.spec_unit.set('a')
        elif selection == 'HARPS(N) S1D NEW DRS':
            self.wavecol.set(1); self.fluxcol.set(3); self.snrcol.set(0); self.errcol.set(4); self.echcol.set(0); self.spec_unit.set('a')
        elif selection == 'HARPS(N) S2D NEW DRS':
            self.wavecol.set(4); self.fluxcol.set(1); self.snrcol.set(0); self.errcol.set(2); self.echcol.set(0); self.spec_unit.set('a')
        else:
            self.wavecol.set(1); self.fluxcol.set(2); self.snrcol.set(0); self.errcol.set(0); self.echcol.set(0); self.spec_unit.set('a')
        self.wavecol_spin.setValue(self.wavecol.get())
        self.fluxcol_spin.setValue(self.fluxcol.get())
        self.snrcol_spin.setValue(self.snrcol.get())
        self.errcol_spin.setValue(self.errcol.get())
        self.echcol_spin.setValue(self.echcol.get())

    def change_units(self, selection):
        self.units_default.set(selection)
        self.spec_unit.set(selection[0])

    def mask_change_units(self, selection):
        self.mask_units_default.set(selection)
        self.mask_unit.set(selection[0])

    def change_frame(self, selection):
        self.wave_frame.set(selection)
        self.spec_vacuum.set(1 if selection == self.wave_values[0] else 0)

    def mask_change_frame(self, selection):
        self.mask_wave_frame.set(selection)
        self.mask_vacuum.set(1 if selection == self.mask_wave_values[0] else 0)

    def change_rvmin(self, selection):
        self.rvmin.set(selection)
        self.rvmin_ext_spin.setValue(selection)
        self.rvmin_prf_spin.setValue(selection)
        
    def change_rvmax(self, selection):
        self.rvmax.set(selection)
        self.rvmax_ext_spin.setValue(selection)
        self.rvmax_prf_spin.setValue(selection)
        
    def change_rvstep(self, selection):
        self.rvstep.set(selection)
        self.rvstep_ext_spin.setValue(selection)
        self.rvstep_prf_spin.setValue(selection)

    def la_change_rvmin(self, selection):
        self.limitlow.set(selection)
        self.limitlow_spin.setValue(selection)
        self.la_limitlow_spin.setValue(selection)
        
    def la_change_rvmax(self, selection):
        self.limitup.set(selection)
        self.limitup_spin.setValue(selection)
        self.la_limitup_spin.setValue(selection)

    def option_la_change(self, selection):
        self.option_limits.set(selection)

    def la_change_fit(self, selection):
        self.option_fit.set(selection)

    def la_change_fit2(self, selection):
        self.option_fit2.set(selection)

    def la_change_fit3(self, selection):
        self.option_fit3.set(selection)
        
    def ts_xdata_change(self, selection):
        self.ts_xplot.set(selection)
        
    def ts_ydata_change(self, selection):
        self.ts_yplot.set(selection)
        
    def ts_comp_change(self, selection):
        self.ts_comp.set(selection)

    def nor_load_indir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select input directory', self.nor_indir.get() or self.basedir)
        if dirname:
            self.nor_indir_edit.setText(dirname)
            self.nor_indir.set(dirname)
            outdir = os.path.join(dirname, self.outdir)
            self.nor_outdir.set(outdir)
            self.nor_outdir_edit.setText(outdir)
            self.prf_indir.set(outdir)
            self.prf_indir_edit.setText(outdir)
            self.prf_outdir.set(outdir)
            self.prf_outdir_edit.setText(outdir)
            self.la_indir.set(outdir)
            self.la_indir_edit.setText(outdir)
            self.la_outdir.set(outdir)
            self.la_outdir_edit.setText(outdir)

    def nor_load_outdir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output directory', self.nor_outdir.get() or self.basedir)
        if dirname:
            self.nor_outdir_edit.setText(dirname)
            self.nor_outdir.set(dirname)
            self.prf_indir_edit.setText(dirname)
            self.prf_indir.set(dirname)
            self.prf_outdir_edit.setText(dirname)
            self.prf_outdir.set(dirname)
            self.la_indir_edit.setText(dirname)
            self.la_indir.set(dirname)
            self.la_outdir_edit.setText(dirname)
            self.la_outdir.set(dirname)

    def nor_load_file(self):
        # load spectrum file
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select spectrum file', self.nor_indir.get() or self.basedir)

        if filename:
            self.nor_spec_edit.setText(filename)
            self.nor_spec.set(filename)
            self.nor_indir.set(os.path.dirname(filename))
            self.nor_indir_edit.setText(os.path.dirname(filename))
            self.nor_outdir.set(os.path.join(self.nor_indir.get(), self.outdir))
            self.nor_outdir_edit.setText(os.path.join(self.nor_indir.get(), self.outdir))
            self.prf_indir.set(self.nor_outdir.get())
            self.prf_indir_edit.setText(self.nor_outdir.get())
            self.prf_outdir.set(self.nor_outdir.get())
            self.prf_outdir_edit.setText(self.nor_outdir.get())
            self.la_indir.set(self.nor_outdir.get())
            self.la_indir_edit.setText(self.nor_outdir.get())
            self.la_outdir.set(self.nor_outdir.get())
            self.la_outdir_edit.setText(self.nor_outdir.get())
        return

    def prf_load_indir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select input directory', self.prf_indir.get() or self.basedir)
        if dirname:
            self.prf_indir_edit.setText(dirname)
            self.prf_indir.set(dirname)
            self.prf_outdir.set(dirname)
            self.prf_outdir_edit.setText(dirname)
            self.la_indir.set(dirname)
            self.la_indir_edit.setText(dirname)
            self.la_outdir.set(dirname)
            self.la_outdir_edit.setText(dirname)

    def prf_load_outdir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output directory', self.prf_outdir.get() or self.basedir)
        if dirname:
            self.prf_outdir_edit.setText(dirname)
            self.prf_outdir.set(dirname)
            self.la_indir_edit.setText(dirname)
            self.la_indir.set(dirname)
            self.la_outdir_edit.setText(dirname)
            self.la_outdir.set(dirname)

    def prf_load_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select normalised file', self.prf_indir.get() or self.basedir)
        if filename:
            self.prf_spec_edit.setText(filename)
            self.prf_spec.set(filename)
            self.prf_indir_edit.setText(os.path.dirname(filename))
            self.prf_indir.set(os.path.dirname(filename))
            

    def mask_load_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select mask file', self.nor_indir.get() or self.basedir)
        if filename:
            self.mask_edit.setText(filename)
            self.mask.set(filename)


    def la_load_indir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select input directory', self.la_indir.get() or self.basedir)
        if dirname:
            self.la_indir_edit.setText(dirname)
            self.la_indir.set(dirname)
            self.la_outdir.set(dirname)
            self.la_outdir_edit.setText(dirname)
            self.ts_indir_edit.setText(dirname)
            self.ts_indir.set(dirname)
            self.ts_outdir_edit.setText(dirname)
            self.ts_outdir.set(dirname)

    def la_load_outdir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output directory', self.la_outdir.get() or self.basedir)
        if dirname:
            self.la_outdir_edit.setText(dirname)
            self.la_outdir.set(dirname)
            self.ts_indir_edit.setText(dirname)
            self.ts_indir.set(dirname)

    def la_load_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select profile file', self.la_indir.get() or self.basedir)
        if filename:
            self.la_spec_edit.setText(filename)
            self.la_spec.set(filename)
            self.la_indir_edit.setText(os.path.dirname(filename))
            self.la_indir.set(os.path.dirname(filename))
            self.ts_spec_edit.setText(filename)
            self.ts_spec.set(filename)
            self.ts_indir_edit.setText(os.path.dirname(filename))
            self.ts_indir.set(os.path.dirname(filename))

    def ts_load_indir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select input directory', self.ts_indir.get() or self.basedir)
        if dirname:
            self.ts_indir_edit.setText(dirname)
            self.ts_indir.set(dirname)
            self.ts_outdir.set(dirname)
            self.ts_outdir_edit.setText(dirname)

    def ts_load_outdir(self):
        dirname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select output directory', self.ts_outdir.get() or self.basedir)
        if dirname:
            self.ts_outdir_edit.setText(dirname)
            self.ts_outdir.set(dirname)

    def ts_load_file(self):
        filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select profile file', self.ts_indir.get() or self.basedir)
        if filename:
            self.ts_spec_edit.setText(filename)
            self.ts_spec.set(filename)
            self.ts_indir_edit.setText(os.path.dirname(filename))
            self.ts_indir.set(os.path.dirname(filename))

    def reset_log(self):
        self.sync_ui_values()
        self.log_text.clear()
        self.log_text.setPlainText(self.messages.get())
        return

    def save_log(self):
        self.sync_ui_values()
        default_name = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S') + '_log_shiva.txt'
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save log file', os.path.join(self.nor_indir.get() or self.basedir, default_name), 'Text Files (*.txt)')
        if fname:
            with open(fname, 'w') as f:
                f.write(self.log_text.toPlainText())

    def thread_normalise(self):
        self.sync_ui_values()
        self.run_in_thread(self.normalise)

    def thread_do_ccf_profile(self):
        self.do_lsd.set(0)
        self.do_ccf.set(1)
        self.sync_ui_values()
        self.run_in_thread(self.do_profile)
        
    def thread_do_lsd_profile(self):
        self.do_lsd.set(1)
        self.do_ccf.set(0)
        self.sync_ui_values()
        self.run_in_thread(self.do_profile)        

    def thread_extract_line(self):
        self.sync_ui_values()
        self.run_in_thread(self.extract_line)

    def thread_norm_profile(self):
        self.sync_ui_values()
        self.run_in_thread(self.norm_profile)

    def thread_fit_profile(self):
        self.sync_ui_values()
        self.run_in_thread(self.fit_profile)

    def thread_moments(self):
        self.sync_ui_values()
        self.run_in_thread(self.moments)

    def thread_bisector(self):
        self.sync_ui_values()
        self.run_in_thread(self.bisector)

    def thread_fourier(self):
        self.sync_ui_values()
        self.run_in_thread(self.fourier)

    def thread_plot_ts(self):
        self.sync_ui_values()
        self.run_in_thread(self.plot_ts)

    def run_in_thread(self, target):
        if self.thread_running.get():
            self.update_log('Another process is already running, wait for it to finish.', timestamp=False)
            return
        self.thread_running.set(1)
        worker = Worker(lambda: self._run_target(target))
        QtCore.QThreadPool.globalInstance().start(worker)

    def _run_target(self, target):
        try:
            target()
        except Exception as e:
            self.update_log(f'Error: {e}', timestamp=False)
        finally:
            self.thread_running.set(0)

    def normalise(self):
        self.update_log('\n#######################', timestamp=False)
        self.update_log('   START normalising the spectra', timestamp=False)
        self.update_log('#######################\n', timestamp=False)
        self.update_log('Parameters:', timestamp=False)
        self.update_log(f'Units: {self.spec_unit.get()}', timestamp=False)
        self.update_log(f'Polynomial degree: {self.degree.get()}', timestamp=False)
        self.update_log(f'Echelle: {self.echcol.get()} (if 0, subsets are used)', timestamp=False)
        self.update_log(f'Subsets: {self.n_ord.get()} (if 0, whole spectrum is used)', timestamp=False)
        self.update_log(f'Refine: {self.refine.get()}\n', timestamp=False)
        spectra = self.find_file(self.nor_spec.get(), self.nor_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.')
            self.update_log('\n#######################', timestamp=False)
            self.update_log('   END normalising the spectra', timestamp=False)
            self.update_log('#######################\n', timestamp=False)
            return
        self.update_log(f'Found {len(spectra)} spectra. Start normalising.\n')
        for n, spec in enumerate(spectra):
            if self.abort_value.get():
                self.update_log("Normalisation aborted.\n")
                self.abort_value.set(False)
                break
            self.update_log(f'{n+1} of {len(spectra)} spectra: working on spectrum {os.path.basename(spec)}.')
            spectrum = pa.read_spectrum(spec, unit=self.spec_unit.get(), wavecol=self.wavecol.get(), fluxcol=self.fluxcol.get(), snrcol=self.snrcol.get(), echcol=self.echcol.get(), errcol=self.errcol.get(), vacuum=bool(self.spec_vacuum.get()))
            nor = pa.norm_spectrum(spectrum['wave'], spectrum['flux'], spectrum['snr'], spectrum['echelle'], deg=self.degree.get(), n_ord=self.n_ord.get(), refine=self.refine.get(), output=False)
            hea = spectrum['header']
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_nor[0], self.degree.get(), self.key_nor[1])
            if self.n_ord.get():
                hea = self.insert_key(hea, self.key_sub[0], self.n_ord.get(), self.key_sub[1])
            else:
                try:
                    subsets = max(spectrum['echelle']) - min(spectrum['echelle']) + 1
                except TypeError:
                    subsets = self.n_ord.get()
                hea = self.insert_key(hea, self.key_sub[0], subsets, self.key_sub[1])
            hea = self.insert_key(hea, self.key_refine[0], bool(self.refine.get()), self.key_refine[1])
            hea = self.insert_key(hea, self.key_swave[0], self.spec_unit.get(), self.key_swave[1])
            self.save_fits(nor, hea, spec, self.nor_outdir.get(), self.nor_output.get())
            self.update_figure(nor['wave'], nor['nflux'], ymin=0, ymax=1.5)
            self.update_log(f'Spectrum {n+1} normalised.')
        self.update_log('\n#######################', timestamp=False)
        self.update_log('   END normalising the spectra', timestamp=False)
        self.update_log('#######################\n', timestamp=False)

    def do_profile(self):
        self.update_log('\n###########################', timestamp=False)
        if self.do_ccf.get():
            self.update_log('  START computing CCF profiles ', timestamp=False)
            root_key = self.base_key['ccf']
        else:
            self.update_log('  START computing LSD profiles ', timestamp=False)
            root_key = self.base_key['lsd']
        self.update_log('###########################\n', timestamp=False)
        self.update_log('Parameters:', timestamp=False)
        self.update_log(f'Spectral units: {self.prf_spec_unit.get()}', timestamp=False)
        self.update_log(f'Absorption mask: {self.mask_invert.get()}', timestamp=False)
        self.update_log(f'Mask units: {self.mask_unit.get()}', timestamp=False)
        self.update_log(f'Mask select line from depths: ({self.mask_dlow.get()},{self.mask_dup.get()})', timestamp=False)
        self.update_log(f'Mask select line from wavelength: ({self.mask_wmin.get()},{self.mask_wmax.get()})', timestamp=False)
        self.update_log(f'Exclude Balmer regions: {bool(self.mask_balmer.get())}', timestamp=False)
        self.update_log(f'Exclude telluric regions: {bool(self.mask_tell.get())}', timestamp=False)
        self.update_log(f'Remove cosmics: {bool(self.cosmic.get())}', timestamp=False)
        self.update_log(f'Clean profile: {bool(self.clean.get())}', timestamp=False)
        self.update_log(f'RV range and step: ({self.rvmin.get()},{self.rvmax.get()}) {self.rvstep.get()}\n', timestamp=False)
        spectra = self.find_file(self.prf_spec.get(), self.prf_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.')
            self.update_log('\n###########################', timestamp=False)
            self.update_log('  END computing profiles ', timestamp=False)
            self.update_log('###########################\n', timestamp=False)
            return
        self.update_log(f'Found {len(spectra)} spectra.')
        try:
            maskfile = self.find_file(self.mask.get())[0]
            #self.update_log(f'{self.mask_cflux.get()}')
            mask = pa.read_mask(maskfile, unit=self.mask_unit.get(), wavecol=int(self.mask_cwave.get()), fluxcol=int(self.mask_cflux.get()), ele=self.mask_els.get(), no_ele=self.mask_noels.get(), depths=(self.mask_dlow.get(), self.mask_dup.get()), balmer=bool(self.mask_balmer.get()), tellurics=bool(self.mask_tell.get()), wmin=self.mask_wmin.get(), wmax=self.mask_wmax.get(), invert=bool(self.mask_invert.get()), vacuum=bool(self.mask_vacuum.get()))
            self.update_log(f'Read mask {os.path.basename(maskfile)}.')
            #self.update_log(f'{mask}.')
            #self.update_figure(mask['wave'], mask['depths'], ymin=None, ymax=None)
            #return
        except Exception:
            self.update_log('No mask was defined. Aborted.')
            self.update_log('\n###########################', timestamp=False)
            self.update_log('  END computing profiles ', timestamp=False)
            self.update_log('###########################\n', timestamp=False)
            return
        for n, spec in enumerate(spectra):
            if self.abort_value.get():
                self.update_log("Profile computation aborted.\n")
                self.abort_value.set(False)
                break
            self.update_log(f'Working on spectrum {n+1} of {len(spectra)}.')
            spectrum = pa.read_spectrum(spec, unit=self.prf_spec_unit.get(), wavecol=self.prf_wavecol.get(), fluxcol=self.prf_fluxcol.get(), snrcol=self.prf_snrcol.get(), nfluxcol=self.prf_nfluxcol.get())
            try:
                for jd_key in self.input_jd:
                    try:
                        jd = spectrum['header'][jd_key]
                        ori_jd_key = jd_key
                        break
                    except KeyError:
                        jd = False
                        ori_jd_key = 'NONE'
            except OSError:
                jd = False
                ori_jd_key = 'NONE'
            
            if self.do_ccf.get():
                profile = pa.compute_ccf(spectrum, mask, vrange=(self.rvmin.get(), self.rvmax.get()), step=self.rvstep.get(), mask_spectrum=bool(self.mask_spectrum.get()), cosmic=bool(self.cosmic.get()), clean=bool(self.clean.get()), weights=bool(self.ccfweight.get()), verbose=False, output=False)
            else:
                profile = pa.compute_lsd(spectrum, mask, vrange=(self.rvmin.get(), self.rvmax.get()), step=self.rvstep.get(), cosmic=bool(self.cosmic.get()), clean=bool(self.clean.get()), verbose=False, output=False)
            hea = fits.PrimaryHDU().header
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_lsdin[0])), os.path.basename(spec), self.key_lsdin[1])
            hea = self.insert_key(hea, self.key_jd[0], jd, ''.join((self.key_jd[1],ori_jd_key)))
            hea = self.insert_key(hea, ' '.join((root_key, self.key_mask[0])), os.path.basename(maskfile), self.key_mask[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_mask_invert[0])), bool(self.mask_invert.get()), self.key_mask_invert[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_mask_spectrum[0])), bool(self.mask_spectrum.get()), self.key_mask_spectrum[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_els[0])), self.mask_els.get(), self.key_els[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_noels[0])), self.mask_noels.get(), self.key_noels[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_dlow[0])), self.mask_dlow.get(), self.key_dlow[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_dup[0])), self.mask_dup.get(), self.key_dup[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_wmin[0])), self.mask_wmin.get(), self.key_wmin[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_wmax[0])), self.mask_wmax.get(), self.key_wmax[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_balmer[0])), bool(self.mask_balmer.get()), self.key_balmer[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_tell[0])), bool(self.mask_tell.get()), self.key_tell[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_cosmic[0])), bool(self.cosmic.get()), self.key_cosmic[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_clean[0])), bool(self.clean.get()), self.key_clean[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_rvmin[0])), self.rvmin.get(), self.key_rvmin[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_rvmax[0])), self.rvmax.get(), self.key_rvmax[1])
            hea = self.insert_key(hea, ' '.join((root_key, self.key_rvstep[0])), self.rvstep.get(), self.key_rvstep[1])
            hea = self.insert_key(hea, self.key_ccfweight[0], bool(self.ccfweight.get()), self.key_ccfweight[1])
            self.save_fits(profile, hea, spec, self.prf_outdir.get(), self.prf_output.get())
            self.update_figure(profile['rv_range'], profile['profile'], ymin=None, ymax=None)
            self.update_log(f'Profile {n+1} computed.')
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END computing profiles ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

    def extract_line(self):
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  START extracting single line ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)
        self.update_log(f'Parameters:\nSpectral units: {self.prf_spec_unit.get()}', timestamp=False)
        self.update_log(f'Central wavelength of line: {self.ext_wave.get()}', timestamp=False)
        self.update_log(f'RV range and step: ({self.rvmin.get()},{self.rvmax.get()}) {self.rvstep.get()}\n', timestamp=False)
        spectra = self.find_file(self.prf_spec.get(), self.prf_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.')
            self.update_log('\n###########################', timestamp=False)
            self.update_log('  END extracting single line ', timestamp=False)
            self.update_log('###########################\n', timestamp=False)
            return
        for n, spec in enumerate(spectra):
            if self.abort_value.get():
                self.update_log("Line extraction aborted.\n")
                self.abort_value.set(False)
                break
            self.update_log(f'Working on spectrum {n+1} of {len(spectra)}.')
            spectrum = pa.read_spectrum(spec, unit=self.prf_spec_unit.get(), wavecol=self.prf_wavecol.get(), fluxcol=self.prf_fluxcol.get(), snrcol=self.prf_snrcol.get(), nfluxcol=self.prf_nfluxcol.get())
            try:
                for jd_key in self.input_jd:
                    try:
                        jd = spectrum['header'][jd_key]
                        ori_jd_key = jd_key
                        break
                    except KeyError:
                        jd = False
                        ori_jd_key = 'NONE'
            except OSError:
                jd = False
                ori_jd_key = 'NONE'
            
            profile = pa.extract_line(spectrum, unit=self.prf_spec_unit.get(), w0=self.ext_wave.get(), vrange=(self.rvmin.get(), self.rvmax.get()), step=self.rvstep.get(), verbose=False, output=False)
            hea = fits.PrimaryHDU().header
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, ' '.join((self.base_key['line'], self.key_lsdin[0])), os.path.basename(spec), self.key_lsdin[1])
            hea = self.insert_key(hea, self.key_jd[0], jd, ''.join((self.key_jd[1],ori_jd_key)))
            hea = self.insert_key(hea, ' '.join((self.base_key['line'], self.key_rvmin[0])), self.rvmin.get(), self.key_rvmin[1])
            hea = self.insert_key(hea, ' '.join((self.base_key['line'], self.key_rvmax[0])), self.rvmax.get(), self.key_rvmax[1])
            hea = self.insert_key(hea, ' '.join((self.base_key['line'], self.key_rvstep[0])), self.rvstep.get(), self.key_rvstep[1])
            self.save_fits(profile, hea, spec, self.prf_outdir.get(), self.ext_output.get())
            self.update_figure(profile['rv_range'], profile['profile'], ymin=None, ymax=None)
            self.update_log(f'Line {n+1} extracted.')
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END extracting single line ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

    def norm_profile(self):
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  START normalising profiles ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)
        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.')
            self.update_log('\n###########################', timestamp=False)
            self.update_log('  END normalising profiles ', timestamp=False)
            self.update_log('###########################\n', timestamp=False)
            return
        if self.std.get():
            stdbasedir = self.la_outdir.get()
            stdbasename = os.path.join(stdbasedir, f"{datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S')}_line_mean_std")
            pfns, stds = pa.norm_profile(spectra, rvcol=1, prfcol=2, errcol=3, sfx=False, std=stdbasename, limits=(self.limitlow.get(), self.limitup.get()))
        else:
            pfns, stds = pa.norm_profile(spectra, rvcol=1, prfcol=2, errcol=3, sfx=False, std=False, limits=(self.limitlow.get(), self.limitup.get()))
        for n, pfn in enumerate(pfns):
            if self.abort_value.get():
                self.update_log("Profile normalisation aborted.\n")
                self.abort_value.set(False)
                break
            hea = pfn.pop('header')
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spectra[n]), self.key_lainp[1])
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lnor[0], True, self.key_lnor[1])
            hea = self.insert_key(hea, self.key_lnrvmin[0], self.limitlow.get(), self.key_lnrvmin[1])
            hea = self.insert_key(hea, self.key_lnrvmax[0], self.limitup.get(), self.key_lnrvmax[1])
            self.save_fits(pfn, hea, spectra[n], self.la_outdir.get(), self.norprf_output.get())
            self.update_figure(pfn['rv_range'], pfn['profile'], y_adds=[pfn.get('nprofile')], ymin=None, ymax=None)
            self.update_log(f'Profile {n+1} normalised.')
        if self.std.get() and stds is not None:
            self.update_figure(stds['rv_mean'], stds['ccf_mean'], y_res=1+stds['std_dev'], ymin=None, ymax=None)
        
        self.ts_spec_edit.setText(f"*{self.norprf_output.get()}")
        self.la_spec_edit.setText(f"*{self.norprf_output.get()}")
        
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END normalising profiles ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

    def fit_profile(self):
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  START fitting profiles ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)
        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.')
            self.update_log('\n###########################', timestamp=False)
            self.update_log('  END fitting profiles ', timestamp=False)
            self.update_log('###########################\n', timestamp=False)
            return
        # Define parameters 
                 
        fit_funs = self.option_fit.get()
        comps = 1
        resolution = self.resolution.get()
        rv0s = str(self.rv0.get())
        widths = str(self.width.get())
        lds = str(self.ld.get())
       
        fit_fun2 = self.option_fit2.get()
        fit_fun3 = self.option_fit3.get()
        
        if fit_fun2 != self.option_la_fits[0]:
            comps = 2
            fit_funs = ','.join((fit_funs,fit_fun2))
            rv0s = ','.join((rv0s, str(self.rv02.get())))
            widths = ','.join((widths, str(self.width2.get())))
            lds = ','.join((lds, str(self.ld2.get())))      
            
            if fit_fun3 != self.option_la_fits[0]:
                comps = 3
                fit_funs = ','.join((fit_funs,fit_fun3))
                rv0s = ','.join((rv0s, str(self.rv03.get())))
                widths = ','.join((widths, str(self.width3.get())))
                lds = ','.join((lds, str(self.ld3.get())))            

        
        for n, spec in enumerate(spectra):
            if self.abort_value.get():
                self.update_log("Profile fitting aborted.")
                self.abort_value.set(False)
                break
            with fits.open(spec) as hdu:
                cols = hdu[1].columns.names
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            flux = data.field(1)
            errs = data.field(2) if self.fit_errs.get() else 0
            
            fit_prf = pa.fit_profile(vrad, flux, errs=errs, fit=fit_funs, rv0=rv0s, width=widths, ld=lds, resolution=resolution, component=comps)

            self.update_log(f'Profile fit results:', timestamp=False)
            
            self.update_log(f"Component 1: Fitting function {fit_prf[0]['function']}", timestamp=False)
            self.update_log(f"RV = {fit_prf[0]['rv']:.4f} +/- {fit_prf[0]['e_rv']:.4f} km/s", timestamp=False)
            self.update_log(f"FWHM/vsini = {fit_prf[0]['width']:.4f} +/- {fit_prf[0]['e_width']:.4f} km/s", timestamp=False)
            self.update_log(f"EW = {fit_prf[0]['ew']:.4f} +/- {fit_prf[0]['e_ew']:.4f} km/s", timestamp=False)
            
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spec), self.key_lainp[1])
            
            hea = self.insert_key(hea, self.key_comps[0], comps, self.key_comps[1])
            hea = self.insert_key(hea, self.key_res[0], resolution, self.key_res[1])
            
            # First component
            hea = self.insert_key(hea, self.key_fit1[0], fit_prf[0]['function'], self.key_fit1[1])
            hea = self.insert_key(hea, self.key_rvguess1[0], self.rv0.get(), self.key_rvguess1[1])
            hea = self.insert_key(hea, self.key_wguess1[0], self.width.get(), self.key_wguess1[1])
            hea = self.insert_key(hea, self.key_fitld1[0], self.ld.get(), self.key_fitld1[1])
            hea = self.insert_key(hea, self.key_rv1[0], fit_prf[0]['rv'], self.key_rv1[1])
            hea = self.insert_key(hea, self.key_rverr1[0], fit_prf[0]['e_rv'], self.key_rverr1[1])
            hea = self.insert_key(hea, self.key_fwhm1[0], fit_prf[0]['width'], self.key_fwhm1[1])
            hea = self.insert_key(hea, self.key_fwhmerr1[0], fit_prf[0]['e_width'], self.key_fwhmerr1[1])
            hea = self.insert_key(hea, self.key_fitew1[0], fit_prf[0]['ew'], self.key_fitew1[1])
            hea = self.insert_key(hea, self.key_fitewerr1[0], fit_prf[0]['e_ew'], self.key_fitewerr1[1])
            
            #Second component
            if fit_fun2 == self.option_la_fits[0]:
                hea = self.insert_key(hea, self.key_fit2[0], 'None', self.key_fit2[1])
                hea = self.insert_key(hea, self.key_rvguess2[0],'None', self.key_rvguess2[1])
                hea = self.insert_key(hea, self.key_wguess2[0], 'None', self.key_wguess2[1])
                hea = self.insert_key(hea, self.key_fitld2[0], 'None', self.key_fitld2[1])
                hea = self.insert_key(hea, self.key_rv2[0], 'None', self.key_rv2[1])
                hea = self.insert_key(hea, self.key_rverr2[0], 'None', self.key_rverr2[1])
                hea = self.insert_key(hea, self.key_fwhm2[0], 'None', self.key_fwhm2[1])
                hea = self.insert_key(hea, self.key_fwhmerr2[0], 'None', self.key_fwhmerr2[1])
                hea = self.insert_key(hea, self.key_fitew2[0], 'None', self.key_fitew2[1])
                hea = self.insert_key(hea, self.key_fitewerr2[0], 'None', self.key_fitewerr2[1]) 
            
            else:
                hea = self.insert_key(hea, self.key_fit2[0], fit_prf[1]['function'], self.key_fit2[1])
                hea = self.insert_key(hea, self.key_rvguess2[0], self.rv02.get(), self.key_rvguess2[1])
                hea = self.insert_key(hea, self.key_wguess2[0], self.width2.get(), self.key_wguess2[1])
                hea = self.insert_key(hea, self.key_fitld2[0], self.ld2.get(), self.key_fitld2[1])
                hea = self.insert_key(hea, self.key_rv2[0], fit_prf[1]['rv'], self.key_rv2[1])
                hea = self.insert_key(hea, self.key_rverr2[0], fit_prf[1]['e_rv'], self.key_rverr2[1])
                hea = self.insert_key(hea, self.key_fwhm2[0], fit_prf[1]['width'], self.key_fwhm2[1])
                hea = self.insert_key(hea, self.key_fwhmerr2[0], fit_prf[1]['e_width'], self.key_fwhmerr2[1])
                hea = self.insert_key(hea, self.key_fitew2[0], fit_prf[1]['ew'], self.key_fitew2[1])
                hea = self.insert_key(hea, self.key_fitewerr2[0], fit_prf[1]['e_ew'], self.key_fitewerr2[1]) 

            #Third component
            if fit_fun3 == self.option_la_fits[0]:
                hea = self.insert_key(hea, self.key_fit3[0], 'None', self.key_fit3[1])
                hea = self.insert_key(hea, self.key_rvguess3[0],'None', self.key_rvguess3[1])
                hea = self.insert_key(hea, self.key_wguess3[0], 'None', self.key_wguess3[1])
                hea = self.insert_key(hea, self.key_fitld3[0], 'None', self.key_fitld3[1])
                hea = self.insert_key(hea, self.key_rv3[0], 'None', self.key_rv3[1])
                hea = self.insert_key(hea, self.key_rverr3[0], 'None', self.key_rverr3[1])
                hea = self.insert_key(hea, self.key_fwhm3[0], 'None', self.key_fwhm3[1])
                hea = self.insert_key(hea, self.key_fwhmerr3[0], 'None', self.key_fwhmerr3[1])
                hea = self.insert_key(hea, self.key_fitew3[0], 'None', self.key_fitew3[1])
                hea = self.insert_key(hea, self.key_fitewerr3[0], 'None', self.key_fitewerr3[1]) 
            
            else:
                hea = self.insert_key(hea, self.key_fit3[0], fit_prf[2]['function'], self.key_fit3[1])
                hea = self.insert_key(hea, self.key_rvguess3[0], self.rv03.get(), self.key_rvguess3[1])
                hea = self.insert_key(hea, self.key_wguess3[0], self.width3.get(), self.key_wguess3[1])
                hea = self.insert_key(hea, self.key_fitld3[0], self.ld3.get(), self.key_fitld3[1])
                hea = self.insert_key(hea, self.key_rv3[0], fit_prf[2]['rv'], self.key_rv3[1])
                hea = self.insert_key(hea, self.key_rverr3[0], fit_prf[2]['e_rv'], self.key_rverr3[1])
                hea = self.insert_key(hea, self.key_fwhm3[0], fit_prf[2]['width'], self.key_fwhm3[1])
                hea = self.insert_key(hea, self.key_fwhmerr3[0], fit_prf[2]['e_width'], self.key_fwhmerr3[1])
                hea = self.insert_key(hea, self.key_fitew3[0], fit_prf[2]['ew'], self.key_fitew3[1])
                hea = self.insert_key(hea, self.key_fitewerr3[0], fit_prf[2]['e_ew'], self.key_fitewerr3[1])             
            

            #res_dict = {'rv_range': vrad, 'profile': flux, 'fit': fit_prf['profile']}
            #self.save_fits(res_dict, hea, spec, self.la_outdir.get(), self.fit_output.get())
            #self.update_fits(spec, hea)
            
            new_data = {}
            for col in cols:
                new_data[col.lower()] = data[col]
            
            new_data['fit'] = fit_prf['profile']
            self.save_fits(new_data, hea, spec, self.la_outdir.get())
            
            if self.save_report.get():
                report_name = os.path.basename(spec)
                report_name = os.path.splitext(report_name)[0]
                report_name = ''.join((report_name,'_fit_report.txt'))
                report_name = os.path.join(self.la_outdir.get(),  report_name)
                with open(report_name, 'w') as fit_report:
                    fit_report.write(fit_prf['report'])
                        
            
            self.update_figure(vrad, flux, y_adds=[fit_prf['profile']], ymin=None, ymax=None)
            self.update_log(f'Fitting done for {os.path.basename(spec)}.\n')
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END fitting profiles ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

    def moments(self):
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  START computing the moments', timestamp=False)
        self.update_log('###########################', timestamp=False)
        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.\n')
            return
        for spec in spectra:
            if self.abort_value.get():
                self.update_log("Moments computation aborted.\n")
                self.abort_value.set(False)
                break
            with fits.open(spec) as hdu:
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            flux = data.field(1)
            errs = data.field(2)
            uselimits = self.option_limits.get()
            if uselimits == self.option_la_values[0]:
                x0 = hea[self.key_rv1[0]]
                fwhm = hea[self.key_fwhm1[0]]
                sigma = fwhm / np.sqrt(8 * np.log(2))
                limits = (x0 - 3 * sigma, x0 + 3 * sigma)
            elif uselimits == self.option_la_values[1]:
                x0 = hea[self.key_rv1[0]]
                vsini = hea[self.key_fwhm1[0]]
                limits = (x0 - vsini, x0 + vsini)
            else:
                #limits = (self.la_limitlow.get(), self.la_limitup.get())
                limits = (self.limitlow.get(), self.limitup.get())
            mom = pa.moments(vrad, flux, errs=errs, limits=limits, normalise=True)
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spec), self.key_lainp[1])
            hea = self.insert_key(hea, self.key_momlim[0], str(limits), self.key_momlim[1])
            
            hea = self.insert_key(hea, self.key_mom0[0], mom['m0'], self.key_mom0[1])
            hea = self.insert_key(hea, self.key_mom0err[0], mom['e_m0'], self.key_mom0err[1])
            hea = self.insert_key(hea, self.key_mom1[0], mom['m1'], self.key_mom1[1])
            hea = self.insert_key(hea, self.key_mom1err[0], mom['e_m1'], self.key_mom1err[1])
            hea = self.insert_key(hea, self.key_mom2[0], mom['m2'], self.key_mom2[1])
            hea = self.insert_key(hea, self.key_mom2err[0], mom['e_m2'], self.key_mom2err[1])
            hea = self.insert_key(hea, self.key_mom3[0], mom['m3'], self.key_mom3[1])
            hea = self.insert_key(hea, self.key_mom3err[0], mom['e_m3'], self.key_mom3err[1])
            hea = self.insert_key(hea, self.key_mom4[0], mom['m4'], self.key_mom4[1])
            hea = self.insert_key(hea, self.key_mom4err[0], mom['e_m4'], self.key_mom4err[1])
            hea = self.insert_key(hea, self.key_momfwhm[0], mom['fwhm'], self.key_momfwhm[1])
            hea = self.insert_key(hea, self.key_momfwhmerr[0], mom['e_fwhm'], self.key_momfwhmerr[1])
            hea = self.insert_key(hea, self.key_momskew[0], mom['skewness'], self.key_momskew[1])
            hea = self.insert_key(hea, self.key_momskewerr[0], mom['e_skewness'], self.key_momskewerr[1])
            hea = self.insert_key(hea, self.key_momkurt[0], mom['kurtosis'], self.key_momkurt[1])
            hea = self.insert_key(hea, self.key_momkurterr[0], mom['e_kurtosis'], self.key_momkurterr[1])
            
            
            #self.save_fits(mom, hea, spec, self.la_outdir.get(), self.mom_output.get())
            
            self.update_fits(spec, hea)
            self.update_figure(vrad, flux, limits=limits)
            
            self.update_log(f'Moments computed for {os.path.basename(spec)}.\n')
            self.update_log(f"m0 (EW) = {round(mom['m0'],4)} +/- {round(mom['e_m0'],4)}  km/s")
            self.update_log(f"m1 (RV) = {round(mom['m1'],4)} +/- {round(mom['e_m1'],4)}  km/s")
            self.update_log(f"m2 (Variance) = {round(mom['m2'],4)} +/- {round(mom['e_m2'],4)}")
            self.update_log(f"m3 (Seed for Skewness) = {round(mom['m3'],4)} +/- {round(mom['e_m3'],4)}")
            self.update_log(f"m4 (Seed for Kurtosis) = {round(mom['m4'],4)} +/- {round(mom['e_m4'],4)}")
            self.update_log(f"FWHM (from m2) = {round(mom['fwhm'],4)} +/- {round(mom['e_fwhm'],4)}")
            self.update_log(f"Skewness (from m3) = {round(mom['skewness'],4)} +/- {round(mom['e_skewness'],4)}")
            self.update_log(f"Kurtosis (from m4) = {round(mom['kurtosis'],4)} +/- {round(mom['e_kurtosis'],4)}")
            
            
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END computing the moments ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

    def bisector(self):
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  START computing the bisector', timestamp=False)
        self.update_log('###########################', timestamp=False)
        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.')
            return
        for spec in spectra:
            if self.abort_value.get():
                self.update_log("Bisector computation aborted.")
                self.abort_value.set(False)
                break
            with fits.open(spec) as hdu:
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            flux = data.field(1)
            uselimits = self.option_limits.get()
            if uselimits == self.option_la_values[0]:
                x0 = hea[self.key_rv1[0]]
                fwhm = hea[self.key_fwhm1[0]]
                sigma = fwhm / np.sqrt(8 * np.log(2))
                limits = (x0 - 3 * sigma, x0 + 3 * sigma)
            elif uselimits == self.option_la_values[1]:
                x0 = hea[self.key_rv1[0]]
                vsini = hea[self.key_fwhm1[0]]
                limits = (x0 - vsini, x0 + vsini)
            else:
                #limits = (self.la_limitlow.get(), self.la_limitup.get())
                limits = (self.limitlow.get(), self.limitup.get())
            bis = pa.bisector(vrad, flux, limits=limits)
            bis_dict = {'rv_range': bis['bisvel'], 'bisector': bis['bisflux'], 'error' : bis['biserr']}
            
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spec), self.key_lainp[1])
            hea = self.insert_key(hea, self.key_bislim[0], str(limits), self.key_bislim[1])
            hea = self.insert_key(hea, self.key_bispan[0], bis['bispan'], self.key_bispan[1])
            hea = self.insert_key(hea, self.key_biserr[0], bis['e_bispan'], self.key_biserr[1])
            
            self.update_fits(spec, hea)

            #self.update_log(f"{bis}")
            self.update_figure(vrad, flux, x_adds=bis_dict['rv_range'], y_adds=[bis_dict['bisector']])
            self.update_log(f'Bisector computed for {os.path.basename(spec)}')
            self.update_log(f"Bisector's span = {round(bis['bispan'],4)} +/- {round(bis['e_bispan'],4)} km/s.")
            
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END computing the bisector', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

    def fourier(self):
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  START computing the Fourier Transform', timestamp=False)
        self.update_log('###########################\n', timestamp=False)
        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.\n')
            return
        for spec in spectra:
            if self.abort_value.get():
                self.update_log("Fourier Transform aborted.\n")
                self.abort_value.set(False)
                break
            with fits.open(spec) as hdu:
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            flux = data.field(1)
            errs = data.field(2) if 'ERR' in data.names else False
            uselimits = self.option_limits.get()
            if uselimits == self.option_la_values[0]:
                x0 = hea[self.key_rv1[0]]
                fwhm = hea[self.key_fwhm1[0]]
                sigma = fwhm / np.sqrt(8 * np.log(2))
                limits = (x0 - 3 * sigma, x0 + 3 * sigma)
            elif uselimits == self.option_la_values[1]:
                x0 = hea[self.key_rv1[0]]
                vsini = hea[self.key_fwhm1[0]]
                limits = (x0 - vsini, x0 + vsini)
            else:
                #limits = (self.la_limitlow.get(), self.la_limitup.get())
                limits = (self.limitlow.get(), self.limitup.get())

            fou = pa.fourier(vrad, flux, errs=errs, limits=limits, ld=self.ld.get())

            fou_dict = {'FFT_fr': fou['FFT_fr'], 'FFT_pow': fou['FFT_pow']}
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spec), self.key_lainp[1])
            hea = self.insert_key(hea, self.key_foulim[0], str(limits), self.key_foulim[1])
            hea = self.insert_key(hea, self.key_fouerr[0], fou['FFT_err'], self.key_fouerr[1])
            hea = self.insert_key(hea, self.key_fouz1[0], fou['zeros'][0], self.key_fouz1[1])
            hea = self.insert_key(hea, self.key_fouz1err[0], fou['e_zeros'][0], self.key_fouz1err[1])
            hea = self.insert_key(hea, self.key_fouz2[0], fou['zeros'][1], self.key_fouz2[1])
            hea = self.insert_key(hea, self.key_fouz2err[0], fou['e_zeros'][1], self.key_fouz2err[1])
            hea = self.insert_key(hea, self.key_fouz3[0], fou['zeros'][2], self.key_fouz3[1])
            hea = self.insert_key(hea, self.key_fouz3err[0], fou['e_zeros'][2], self.key_fouz3err[1])
            hea = self.insert_key(hea, self.key_fouratio[0], fou['ratio'], self.key_fouratio[1])
            hea = self.insert_key(hea, self.key_fouratioerr[0], fou['e_ratio'], self.key_fouratioerr[1])
            hea = self.insert_key(hea, self.key_fouvsini1[0], fou['vsini'][0], self.key_fouvsini1[1])
            hea = self.insert_key(hea, self.key_fouvsini1err[0], fou['e_vsini'][0], self.key_fouvsini1err[1])
            hea = self.insert_key(hea, self.key_fouvsini2[0], fou['vsini'][1], self.key_fouvsini2[1])
            hea = self.insert_key(hea, self.key_fouvsini2err[0], fou['e_vsini'][1], self.key_fouvsini2err[1])
            hea = self.insert_key(hea, self.key_fouvsini3[0], fou['vsini'][2], self.key_fouvsini3[1])
            hea = self.insert_key(hea, self.key_fouvsini3err[0], fou['e_vsini'][2], self.key_fouvsini3err[1])
            hea = self.insert_key(hea, self.key_fouvsini[0], fou['mean_vsini'], self.key_fouvsini[1])
            hea = self.insert_key(hea, self.key_fouvsinierr[0], fou['e_mean_vsini'], self.key_fouvsinierr[1])
            hea = self.insert_key(hea, self.key_fourv[0], fou['rv'], self.key_fourv[1])
            
            self.update_fits(spec, hea)
            self.save_fits(fou_dict, hea, spec, self.la_outdir.get(), self.fou_output.get())
                        
            y_add = np.ones(fou_dict['FFT_fr'].shape)*fou['FFT_err']
            if len(fou['zeros']) > 3:
                xmax = fou['zeros'][3]
                max_idx = np.searchsorted(fou_dict['FFT_fr'], xmax)
                self.update_figure(fou_dict['FFT_fr'][:max_idx], fou_dict['FFT_pow'][:max_idx], y_adds=[y_add[:max_idx]], limits=fou['zeros'][:3], logscale=True)
            else:
                self.update_figure(fou_dict['FFT_fr'], fou_dict['FFT_pow'], y_adds=[y_add], limits=fou['zeros'], logscale=True)
            
            self.update_log(f'Fourier transform computed for {os.path.basename(spec)}')
            self.update_log(f"First zero vsini = {np.round(fou['vsini'][0],4)} +/- {np.round(fou['e_vsini'][0],4)} km/s", timestamp=False)
            self.update_log(f"Mean vsini = {np.round(fou['mean_vsini'],4)} +/- {np.round(fou['e_mean_vsini'],4)} km/s", timestamp=False)
            self.update_log(f"q2/q1 = {np.round(fou['ratio'],4)} +/- {np.round(fou['e_ratio'],4)} km/s", timestamp=False)
            self.update_log(f"RV = {np.round(fou['rv'],4)} km/s", timestamp=False)
        
        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END computing the Fourier Transform ', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

    def plot_ts(self):
        #self.option_ts_plot = ['RV', 'EW', 'width', 'm0 (EW)', 'm1 (RV)', 'm2 (sigma)', 'skewness', 'kurtosis', 'bispan', 'vsini (Fourier)', 'q2/q1 (Fourier)']
        #self.option_ts_comp = ['1', '2', '3']

        self.update_log('\n###########################', timestamp=False)
        self.update_log('  START plotting time series', timestamp=False)
        self.update_log('###########################\n', timestamp=False)
        
        xdata = self.ts_xplot.get()
        ydata = self.ts_yplot.get()
        error = bool(self.ts_error.get())
        comp = self.ts_comp.get()
        save = bool(self.ts_save.get())

        self.update_log(f'Plotting {xdata} vs. {ydata} - component {comp}')
        self.update_log(f'With errors: {error}', timestamp=False)

        kwds = {f'{self.option_ts_plot[0]}1' : (self.key_jd[0], None), \
                f'{self.option_ts_plot[0]}2' : (self.key_jd[0], None), \
                f'{self.option_ts_plot[0]}3' : (self.key_jd[0], None), \
                f'{self.option_ts_plot[1]}1' : (self.key_rv1[0], self.key_rverr1[0]), \
                f'{self.option_ts_plot[1]}2' : (self.key_rv2[0], self.key_rverr2[0]), \
                f'{self.option_ts_plot[1]}3' : (self.key_rv3[0], self.key_rverr3[0]), \
                f'{self.option_ts_plot[2]}1' : (self.key_fitew1[0], self.key_fitewerr1[0]), \
                f'{self.option_ts_plot[2]}2' : (self.key_fitew2[0], self.key_fitewerr2[0]), \
                f'{self.option_ts_plot[2]}3' : (self.key_fitew3[0], self.key_fitewerr3[0]), \
                f'{self.option_ts_plot[3]}1' : (self.key_fwhm1[0], self.key_fwhmerr1[0]), \
                f'{self.option_ts_plot[3]}2' : (self.key_fwhm2[0], self.key_fwhmerr2[0]), \
                f'{self.option_ts_plot[3]}3' : (self.key_fwhm3[0], self.key_fwhmerr3[0]), \
                f'{self.option_ts_plot[4]}1' : (self.key_mom0[0], self.key_mom0err[0]), \
                f'{self.option_ts_plot[4]}2' : (self.key_mom0[0], self.key_mom0err[0]), \
                f'{self.option_ts_plot[4]}3' : (self.key_mom0[0], self.key_mom0err[0]), \
                f'{self.option_ts_plot[5]}1' : (self.key_mom1[0], self.key_mom1err[0]), \
                f'{self.option_ts_plot[5]}2' : (self.key_mom1[0], self.key_mom1err[0]), \
                f'{self.option_ts_plot[5]}3' : (self.key_mom1[0], self.key_mom1err[0]), \
                f'{self.option_ts_plot[6]}1' : (self.key_mom2[0], self.key_mom2err[0]), \
                f'{self.option_ts_plot[6]}2' : (self.key_mom2[0], self.key_mom2err[0]), \
                f'{self.option_ts_plot[6]}3' : (self.key_mom2[0], self.key_mom2err[0]), \
                f'{self.option_ts_plot[7]}1' : (self.key_momskew[0], self.key_momskewerr[0]), \
                f'{self.option_ts_plot[7]}2' : (self.key_momskew[0], self.key_momskewerr[0]), \
                f'{self.option_ts_plot[7]}3' : (self.key_momskew[0], self.key_momskewerr[0]), \
                f'{self.option_ts_plot[8]}1' : (self.key_momkurt[0], self.key_momkurterr[0]), \
                f'{self.option_ts_plot[8]}2' : (self.key_momkurt[0], self.key_momkurterr[0]), \
                f'{self.option_ts_plot[8]}3' : (self.key_momkurt[0], self.key_momkurterr[0]), \
                f'{self.option_ts_plot[9]}1' : (self.key_bispan[0], self.key_biserr[0]), \
                f'{self.option_ts_plot[9]}2' : (self.key_bispan[0], self.key_biserr[0]), \
                f'{self.option_ts_plot[9]}3' : (self.key_bispan[0], self.key_biserr[0]), \
                f'{self.option_ts_plot[10]}1' : (self.key_fouvsini1[0], self.key_fouvsini1err[0]), \
                f'{self.option_ts_plot[10]}2' : (self.key_fouvsini1[0], self.key_fouvsini1err[0]), \
                f'{self.option_ts_plot[10]}3' : (self.key_fouvsini1[0], self.key_fouvsini1err[0]), \
                f'{self.option_ts_plot[11]}1' : (self.key_fouratio[0], self.key_fouratioerr[0]), \
                f'{self.option_ts_plot[11]}2' : (self.key_fouratio[0], self.key_fouratioerr[0]), \
                f'{self.option_ts_plot[11]}3' : (self.key_fouratio[0], self.key_fouratioerr[0]) \
                }
                
        key_xdata = kwds[f'{xdata}{comp}']
        key_ydata = kwds[f'{ydata}{comp}']
        #self.update_log(f"{key_xdata}", timestamp=False)
        #self.update_log(f"{key_ydata}", timestamp=False)
        
        spectra = self.find_file(self.ts_spec.get(), self.ts_indir.get())
        if not spectra:
            self.update_log('No input found. Aborted.\n')
            return
        xvalues = np.arange(len(spectra), dtype='float')
        yvalues = np.arange(len(spectra), dtype='float')
        xerrvalues = np.zeros(len(spectra), dtype='float')
        yerrvalues = np.zeros(len(spectra), dtype='float')
        
        print_xlog = False
        print_ylog = False
            
        for n,spec in enumerate(spectra):
            if self.abort_value.get():
                self.update_log("Plotting time series aborted.\n")
                self.abort_value.set(False)
                break
            with fits.open(spec) as hdu:
                hea = hdu[0].header
            if hea[key_xdata[0]]:
                xvalues[n] = hea[key_xdata[0]]
            else:
                print_xlog = True

            if hea[key_ydata[0]]:
                yvalues[n] = hea[key_ydata[0]]
            else:
                print_ylog = True
            if error:
                try:
                    xerrvalues[n] = hea[key_xdata[1]]
                except ValueError:
                    pass
                try:
                    yerrvalues[n] = hea[key_ydata[1]]
                except ValueError:
                    pass
        if print_xlog:
            self.update_log(f"No values found for {xdata}, incremental numbers will be used.", timestamp=False)
        if print_ylog:
            self.update_log(f"No values found for {ydata}, incremental numbers will be used.", timestamp=False)
        
        self.update_figure(xvalues, yvalues, x_err=xerrvalues, y_err=yerrvalues, ymin=None, ymax=None)
        
        if save:
            savedata = np.vstack((xvalues, xerrvalues, yvalues, yerrvalues))            
            outname = f"{datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S')}_{xdata.replace('/','')}_vs_{ydata.replace('/','')}.txt"
            outname = os.path.join(self.ts_outdir.get(), outname)
            header = f'#1. {xdata} - 2. err. {xdata} - 3. {ydata} - 4. err. {ydata}'
            np.savetxt(outname, np.transpose(savedata), header=header)

        self.update_log('\n###########################', timestamp=False)
        self.update_log('  END plotting time series', timestamp=False)
        self.update_log('###########################\n', timestamp=False)

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    app.setStyleSheet(f'QGroupBox {{font: bold;}}')
    window = ShivaQtApp()
    window.show()
    app.exec()
