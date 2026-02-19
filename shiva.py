"""
SHIVA: a Simple and Helpful Interface for Variability Analysis
SHIVA is a GUI wrapper for PARVATI: Profiles Analysis and Radial Velocities using Astronomical Tools for Investigation (a Python package to compute and analyse stellar mean line profiles)
Written by Monica Rainer

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

import os
import time
import datetime

import glob
import warnings

#import sys
#sys.path.append("/home/monica/Documents/GitHub/shiva")
#import test_parvati as pa
import parvati as pa

import numpy as np
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', AstropyWarning)

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)

from threading import Thread

import tkinter.filedialog as filedialog
import tkinter as tk
from tkinter import ttk

from tktooltip import ToolTip

__version__ = "1.0.2"


############
#  GUI     #
############


class MainApp(ttk.Frame):
    def __init__(self, master):
        self.master = master
        self.basedir = os.getcwd()  # root directory of the script
        #self.basedir = '/home/monica/Documents/GAPS/YO02/s1d'
        self.version = ''.join(('SHIVA v',__version__))
        try: self.pa_version = ''.join(('PARVATI v',pa.__version__))
        except AttributeError:
            self.pa_version = ''.join(('PARVATI v','XXX'))
        self.key_prg = 'SP' #SHIVA-PARVATI
        self.sp_logo = 'sp_logo.png'
        self.outdir = 'shiva_output'
        self.font_type = 'Calibri'
        self.label_size = 12
        self.text_size = 10
        self.centreWindow()
        ttk.Frame.__init__(self, self.master)
        self.master.title('SHIVA - Simple and Helpful Interface for Variability Analysis using PARVATI')
        self.abort_value = False
        self.define_tooltip()
        self.define_keys()
        self.define_entries()
        self.set_entries()
        self.define_plot()
        self.create_widgets()


    def centreWindow(self):
        # Define GUI window size and position
        # normalisation widget
        self.norrows = 6
        self.norcols = 17
        # global profile widget
        self.lsdrows = 11
        # profile input sub-widget
        self.lsdrows_input = 3
        # extraction profile sub-widget
        self.extrows = 2
        # creation profile sub-widget
        self.prfrows = 6
        # plot widget
        self.pltrows = 8
        self.pltcols = 11
        # message widget
        self.msgrows = 9
        # global line analysis widget
        self.larows = 14
        self.lacols = 17
        # line analysis input sub-widget
        self.larows_input = 3
        # profile normalisation sub-widget
        self.norprfrows = 2
        # profile fitting sub-widget
        self.fitrows = 3
        # moments sub-widget
        self.momrows = 2
        # bisector sub-widget
        self.bisrows = 2
        # Fourier transform sub-widget
        self.fourows = 2
        # bottom sub-widget
        self.bottomrows = 3
        # whole GUI
        self.all_rows=self.norrows+self.lsdrows
        self.all_columns=self.norcols+self.pltcols+self.lacols

        w = 1600
        h = 850
        sw = self.master.winfo_screenwidth()
        w = int(min(w, sw-100))
        sh = self.master.winfo_screenheight()
        h = int(min(h, sh-100))
        x = int((sw-w)/2)
        y = int((sh-h)/2)
        self.master.geometry('%dx%d+%d+%d' % (w,h,x,y))
        #self.master.geometry('500x300')
        for n in range(self.all_rows):
            self.master.rowconfigure(n, weight=1, minsize=1)
        for n in range(self.all_columns):
            self.master.columnconfigure(n, weight=1, minsize=1)
        #for n in range(self.norcols):
        #    self.master.columnconfigure(n, weight=3, minsize=3)
        #for n in range(self.norcols,self.all_columns):
        #    self.master.columnconfigure(n, weight=1, minsize=3)
        #self.master.columnconfigure(0, weight=6, minsize=3)
        #self.master.columnconfigure(1, weight=1, minsize=3)
        #self.master.columnconfigure(2, weight=4, minsize=3)
        #self.grid_columnconfigure(0,  weight=1)
        #self.grid_rowconfigure(0,  weight=1)

    def define_tooltip(self):
        self.t_delay=0.5
        self.t_tcolor='#000000'
        self.t_background="#ffffe0"
        self.t_borderwidth=2
        self.t_relief='raised'
        self.t_show_duration=5
    
    def define_keys(self):
        
        # Possible JD keywords in input FITS spectra
        self.input_jd = ['MJD-OBS']
        
        # shiva version
        self.key_ver = (f'HIERARCH SHIVA VERSION', 'SHIVA version used')
        # parvati version
        self.pa_ver = (f'HIERARCH PARVATI VERSION', 'PARVATI version used')
        
        # normalisation
        self.key_nor = (f'HIERARCH {self.key_prg} NOR DEG', 'Polynomial degree for normalisation.')
        self.key_sub = (f'HIERARCH {self.key_prg} NOR SUBSETS', 'Subsets independently normalised.')
        self.key_refine = (f'HIERARCH {self.key_prg} NOR REFINE', 'Refined normalisation')
        
        # LSD/CCF
        self.key_jd = (f'HIERARCH {self.key_prg} JD', 'JD value in inputs spectrum ')
        self.key_swave = (f'HIERARCH {self.key_prg} WAVE UNITS', 'Wavelength unit (a=Angstrom, n=nanometer, m=micron)')
        # HIERARCH {self.key_prg} LSD  or HIERARCH {self.key_prg} CCF
        self.base_key = {'lsd' : f'HIERARCH {self.key_prg} LSD', 'ccf' : f'HIERARCH {self.key_prg} CCF', 'line' : f'HIERARCH {self.key_prg} EXT'}
        self.key_lsdin = ('INPUT', 'Input spectrum used to compute profile')
        self.key_mask = ('MASK', 'Mask used to compute profile')
        self.key_mask_invert = ('MASK INVERT', 'Absorption mask, inverted before being used')
        self.key_mask_spectrum = ('MASK SPECTRUM', 'Spectrum (observed or template) used as mask - with continuum')
        # mask invert? mask spectrum?
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
        self.key_ccfweight = (' '.join((self.base_key['ccf'],'WEIGHTED')), 'Weighted CCF, using the S/N')
        
        # Generic Line Analysis
        self.key_lainp = (f'HIERARCH {self.key_prg} LA INPUT', 'Input file for the line analysis')
        self.key_lnor = (f'HIERARCH {self.key_prg} PRF NOR', 'Normalised profile')
        
        # Normalise profile / Manual profile limits
        self.key_lnrvmin = (f'HIERARCH {self.key_prg} LSD NOR RVMIN', 'Minimum RV for continuum definition (km/s)')
        self.key_lnrvmax = (f'HIERARCH {self.key_prg} LSD NOR RVMAX', 'Maximum RV for continuum definition (km/s)')
        
        # Fit the profile
        self.key_fitgauss = (f'HIERARCH {self.key_prg} FIT GAUSS', 'Gaussian fit performed')
        self.key_fitrot = (f'HIERARCH {self.key_prg} FIT ROT', 'Rotational fit performed')
        self.key_fitld = (f'HIERARCH {self.key_prg} FIT LD', 'Linear Limb Darkening')
        self.key_rvguess = (f'HIERARCH {self.key_prg} FIT RV GUESS', 'RV guess value')
        self.key_wguess = (f'HIERARCH {self.key_prg} FIT WIDTH GUESS', 'Line width guess value')
        self.key_gaussrv = (f'HIERARCH {self.key_prg} FIT GAUSS RV', 'RV Gaussian fit')
        self.key_gaussrverr = (f'HIERARCH {self.key_prg} FIT GAUSS RVERR', 'RV error Gaussian fit')
        self.key_gaussfwhm = (f'HIERARCH {self.key_prg} FIT GAUSS FWHM', 'FWHM Gaussian fit')
        self.key_gaussfwhmerr = (f'HIERARCH {self.key_prg} FIT GAUSS FWHMERR', 'FWHM error Gaussian fit')
        self.key_gaussew = (f'HIERARCH {self.key_prg} FIT GAUSS EW', 'EW Gaussian fit')
        self.key_gaussewerr = (f'HIERARCH {self.key_prg} FIT GAUSS EWERR', 'EW error Gaussian fit')
        self.key_rotrv = (f'HIERARCH {self.key_prg} FIT ROT RV', 'RV Rotational fit')
        self.key_rotrverr = (f'HIERARCH {self.key_prg} FIT ROT RVERR', 'RV error Rotational fit')
        self.key_rotvsini = (f'HIERARCH {self.key_prg} FIT ROT VSINI', 'Vsini Rotational fit')
        self.key_rotvsinierr = (f'HIERARCH {self.key_prg} FIT ROT VSINIERR', 'Vsini error Rotational fit')
        self.key_rotew = (f'HIERARCH {self.key_prg} FIT ROT EW', 'EW Rotational fit')
        self.key_rotewerr = (f'HIERARCH {self.key_prg} FIT ROT EWERR', 'EW error Rotational fit')        
        self.key_lorentzrv = (f'HIERARCH {self.key_prg} FIT LORENTZ RV', 'RV Lorentzian fit')
        self.key_lorentzrverr = (f'HIERARCH {self.key_prg} FIT LORENTZ RVERR', 'RV error Lorentzian fit')
        self.key_lorentzfwhm = (f'HIERARCH {self.key_prg} FIT LORENTZ FWHM', 'FWHM Lorentzian fit')
        self.key_lorentzfwhmerr = (f'HIERARCH {self.key_prg} FIT LORENTZ FWHMERR', 'FWHM error Lorentzian fit')
        self.key_lorentzew = (f'HIERARCH {self.key_prg} FIT LORENTZ EW', 'EW Lorentzian fit')
        self.key_lorentzewerr = (f'HIERARCH {self.key_prg} FIT LORENTZ EWERR', 'EW error Lorentzian fit')
        self.key_voigtrv = (f'HIERARCH {self.key_prg} FIT VOIGT RV', 'RV Voigt fit')
        self.key_voigtrverr = (f'HIERARCH {self.key_prg} FIT VOIGT RVERR', 'RV error Voigt fit')
        self.key_voigtfwhm = (f'HIERARCH {self.key_prg} FIT VOIGT FWHM', 'FWHM Voigt fit')
        self.key_voigtfwhmerr = (f'HIERARCH {self.key_prg} FIT VOIGT FWHMERR', 'FWHM error Voigt fit')
        self.key_voigtew = (f'HIERARCH {self.key_prg} FIT VOIGT EW', 'EW Voigt fit')
        self.key_voigtewerr = (f'HIERARCH {self.key_prg} FIT VOIGT EWERR', 'EW error Voigt fit')
        
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
        # Define empty tabs and input variable/buttons
        
        # Normalise spectrum
        self.nor_indir = tk.StringVar()     # input folder
        self.nor_spec = tk.StringVar()      # input spectrum OR pattern
        self.option_instr = tk.StringVar()  # choose the spectrograph
        self.units_default = tk.StringVar()  # wavelength unit (a/n/m)
        self.wave_frame = tk.StringVar()    # wave reference frame (vacuum/air)
        self.wavecol = tk.IntVar()          # ASCII column/FITS table field of the wavelength
        self.fluxcol = tk.IntVar()          # ASCII column/FITS table field of the flux
        self.snrcol = tk.IntVar()           # ASCII column/FITS table field of the SNR (if any)
        self.errcol = tk.IntVar()           # ASCII column/FITS table field of the errors (if any)
        self.echcol = tk.IntVar()           # ASCII column/FITS table field of the echelle orders (if any)
        self.degree = tk.IntVar()           # degree of the normalisation polynomial
        self.n_ord = tk.IntVar()            # number of "fake" echelle orders (if echcol is not defined)
        self.refine = tk.IntVar()           # do we try to refine the normalisation?
        self.norplot= tk.IntVar()           # plot the output
        self.spec_unit = tk.StringVar()     # input spectrum wavelength unit
        self.spec_vacuum = tk.IntVar()      # input spectrum wavelength is in vacuum
        self.nor_outdir = tk.StringVar()    # folder of the output normalised spectra (if missing, no output saved)
        self.nor_output = tk.StringVar()    # suffix of the output normalised spectra (if missing, no output saved)
        
        # Compute profile - INPUT
        self.lsd_indir = tk.StringVar()     # input folder
        self.lsd_spec = tk.StringVar()      # input spectrum OR pattern
        self.lsd_wavecol = tk.IntVar()      # ASCII column/FITS table field of the wavelength
        self.lsd_fluxcol = tk.IntVar()      # ASCII column/FITS table field of the flux
        self.lsd_nfluxcol = tk.IntVar()     # ASCII column/FITS table field of the normalised flux
        self.lsd_snrcol = tk.IntVar()       # ASCII column/FITS table field of the SNR (if any)
        self.lsd_spec_unit = tk.StringVar() # input spectrum wavelength unit
        self.lsd_outdir = tk.StringVar()    # folder of the output LSD profiles (if missing, no output saved)
        # Compute profile - extraction
        self.ext_wave = tk.DoubleVar()      # wavelength of single line
        self.rvmin = tk.DoubleVar()         # minimum RV
        self.rvmax = tk.DoubleVar()         # maximum RV
        self.rvstep = tk.DoubleVar()        # RV step
        self.ext_output = tk.StringVar()    # suffix of the output extracted profiles (if missing, no output saved)
        self.extplot = tk.IntVar()          # plot the output
        # Compute profile - LSD/CCF
        self.mask = tk.StringVar()          # input mask
        self.mask_invert = tk.IntVar()      # invert the mask (from flux to depths)
        self.mask_spectrum = tk.IntVar()    # use a spectrum as mask (with continuum)
        self.mask_unit = tk.StringVar()     # input mask wavelength unit
        self.mask_units_default = tk.StringVar() # input mask wavelength unit (default value)
        self.mask_wave_frame = tk.StringVar()  # mask wave frame combobox
        self.mask_vacuum = tk.IntVar()      # input mask wavelength is in vacuum
        self.mask_dlow = tk.DoubleVar()     # minimum line depths
        self.mask_dup = tk.DoubleVar()      # maximum line depths
        self.mask_wmin = tk.DoubleVar()     # minimum mask wavelength
        self.mask_wmax = tk.DoubleVar()     # maximum mask wavelength
        self.mask_els = tk.StringVar()      # input mask elements to be used
        self.mask_noels = tk.StringVar()    # input mask elements to be excluded
        self.mask_balmer = tk.IntVar()      # exclude Balmer regions
        self.mask_tell = tk.IntVar()        # exclude telluric regions
        self.cosmic = tk.IntVar()           # remove cosmic rays from the spectrum
        self.clean = tk.IntVar()            # clean the spectrum (smoothing spline)
        self.ccfweight = tk.IntVar()        # Weighted CCF using the S/N field
        self.do_lsd = tk.IntVar()           # Compute LSD profiles
        self.do_ccf = tk.IntVar()           # Compute CCF profiles
        self.lsdplot = tk.IntVar()          # plot the output
        self.lsd_output = tk.StringVar()    # suffix of the output LSD profiles (if missing, no output saved)
        
        # MESSAGES
        self.messages = tk.StringVar()      # messages
        
        # Line Analysis
        self.la_indir = tk.StringVar()      # input folder
        self.la_spec = tk.StringVar()       # input (normalised) LSD profile OR pattern
        self.la_outdir = tk.StringVar()     # output folder
        self.option_limits = tk.StringVar()  # choose the line limits definition
        self.show_limitlow  = tk.StringVar()    # show either manual limit or 3sigma/vsini
        self.show_limitup  = tk.StringVar()     # show either manual limit or 3sigma/vsini
        
        # Normalise LSD profiles
        self.limitlow = tk.DoubleVar()      # lower line limit (default self.rvmin + something)
        self.limitup = tk.DoubleVar()       # upper line limit (default self.rvmax - something)
        self.std = tk.IntVar()              # create plot with mean profile and stdev
        self.norprfplot= tk.IntVar()        # plot the output
        self.norprf_output = tk.StringVar() # suffix of the output normalised LSD profiles (if missing, no output saved)
        
        # Fit LSD profiles
        self.fitgauss = tk.IntVar()         # Gaussian fit
        self.fitrot = tk.IntVar()           # rotational fit
        self.fitlorentz = tk.IntVar()       # Lorentzian fit
        self.fitvoigt = tk.IntVar()         # Voigt fit
        self.rv0 = tk.DoubleVar()           # guess RV value
        self.width = tk.DoubleVar()         # guess line width value
        self.ld = tk.DoubleVar()            # linear limb darkening value
        self.fit_errs = tk.IntVar()         # use the data errors when fitting
        self.fitplot= tk.IntVar()           # plot the output
        self.fit_output = tk.StringVar()    # suffix of the output fitting LSD profiles (if missing, no output saved)
        
        # Compute moments
        self.usegauss = tk.IntVar()         # use Gaussian 3sigma to define the line limit
        self.userot = tk.IntVar()           # use rotational width to define the line limit
        self.momplot= tk.IntVar()           # plot the output
        self.mom_output = tk.StringVar()    # suffix of the moments output (if missing, no output saved)
        
        # Compute bisector
        self.bisplot= tk.IntVar()           # plot the output
        self.bis_output = tk.StringVar()    # suffix of the bisector output (if missing, no output saved)
                
        # Compute Fourier transform
        self.fouplot= tk.IntVar()           # plot the output
        self.fou_output = tk.StringVar()    # suffix of the Fourier output (if missing, no output saved)
        
        # Threads
        self.thread_running = tk.IntVar()    # check if threads are running

    def set_entries(self):

        # Define default values for functions

        # COMBOBOX
        # norframe
        self.option_values = ['UNDEF', 'ESPRESSO S1D', 'ESPRESSO S2D', 'GIANO-B MS1D', 'CARMENES', 'HARPS(N) S1D NEW DRS', 'HARPS(N) S2D NEW DRS'] # list of options for the instruments list -- IMPORTANT: if the list is changed, then change accordingly the function change_cols
        self.units_values = ['angstroms','nm','micron']
        self.wave_values = ['Vacuum', 'Air']
        #prfframe
        self.mask_units_values = ['angstroms','nm','micron']
        self.mask_wave_values = ['Vacuum', 'Air']
        self.option_la_values = ['Limits: Gaussian', 'Limits: Rotational', 'Limits: Manual']
        
        # Normalise spectrum
        self.nor_indir.set(self.basedir)    # input folder
        self.nor_spec.set('*.fits')         # input spectrum OR pattern
        self.option_instr.set(self.option_values[0])  # choose the spectrograph
        self.units_default.set(self.units_values[0])  # wavelength unit (a/n/m)
        self.wave_frame.set(self.wave_values[0])    # wave reference frame (vacuum/air)
        self.wavecol.set(1)                 # ASCII column/FITS table field of the wavelength
        self.fluxcol.set(2)                 # ASCII column/FITS table field of the flux
        self.snrcol.set(0)                  # ASCII column/FITS table field of the SNR (if any)
        self.errcol.set(0)                  # ASCII column/FITS table field of the error (if any)
        self.echcol.set(0)                  # ASCII column/FITS table field of the echelle orders (if any)
        self.degree.set(2)                  # degree of the normalisation polynomial
        self.n_ord.set(0)                   # number of "fake" echelle orders (if echcol is not defined)
        self.refine.set(0)                  # do we try to refine the normalisation?
        self.norplot.set(1)                 # plot the output
        self.spec_unit.set('a')             # input spectrum wavelength unit
        self.spec_vacuum.set(1)             # input spectrum wavelength is in vacuum
        self.nor_outdir.set(os.path.join(self.basedir,self.outdir))    # folder of the output normalised spectra (if missing, no output saved)
        self.nor_output.set('_nor.fits')    # suffix of the output normalised spectra (if missing, no output saved)
        
        # Compute profile - INPUT
        self.lsd_indir.set(self.nor_outdir.get())    # input folder
        self.lsd_spec.set('*_nor.fits')     # input normalised spectrum OR pattern
        self.lsd_wavecol.set(1)             # ASCII column/FITS table field of the wavelength
        self.lsd_fluxcol.set(2)             # ASCII column/FITS table field of the flux
        self.lsd_nfluxcol.set(3)            # ASCII column/FITS table field of the normalised flux
        self.lsd_snrcol.set(4)              # ASCII column/FITS table field of the SNR (if any)
        self.lsd_spec_unit.set('a')         # input spectrum wavelength unit
        self.lsd_outdir.set(self.nor_outdir.get())    # folder of the output LSD profiles (if missing, no output saved)
        # Compute profile - extraction
        self.ext_wave.set(6562.801)         # wavelength of single line
        self.rvmin.set(-100)                # minimum RV
        self.rvmax.set(100)                 # maximum RV
        self.rvstep.set(1)                  # RV step
        self.ext_output.set('_ext.fits')    # suffix of the output extracted profiles (if missing, no output saved)
        self.extplot.set(1)                 # plot the output
        # Compute profile - LSD/CCF
        self.mask.set('')                   # input mask
        self.mask_invert.set(0)             # invert the mask (from flux to depths)
        self.mask_spectrum.set(0)           # use a spectrum as mask (with continuum)
        self.mask_unit.set('a')             # input mask wavelength unit
        self.mask_units_default.set(self.mask_units_values[0])
        self.mask_wave_frame.set(self.wave_values[0])  # mask wavelength frame (vacuum/air)      
        self.mask_vacuum.set(1)             # input mask wavelength is in vacuum
        self.mask_dlow.set(0.01)            # minimum line depths
        self.mask_dup.set(1)                # maximum line depths
        self.mask_wmin.set(0)               # minimum mask wavelength
        self.mask_wmax.set(0)               # maximum mask wavelength
        self.mask_els.set('')               # input mask elements to be used
        self.mask_noels.set('')             # input mask elements to be excluded
        self.mask_balmer.set(1)             # exclude Balmer regions
        self.mask_tell.set(1)               # exclude telluric regions
        self.cosmic.set(0)                  # remove cosmic rays from the spectrum
        self.clean.set(0)                   # clean the spectrum (smoothing spline)
        self.ccfweight.set(0)               # Weighted CCF using the S/N field
        self.do_lsd.set(0)                  # Compute LSD profiles
        self.do_ccf.set(1)                  # Compute CCF profiles
        self.lsdplot.set(1)                 # plot the output
        if self.do_lsd.get():
            self.lsd_output.set('_lsd.fits')    # suffix of the output LSD profiles (if missing, no output saved)
        else:
            self.lsd_output.set('_ccf.fits')    # suffix of the output LSD profiles (if missing, no output saved)
            
        # MESSAGES
        self.messages.set(f"{self.version}\n{self.pa_version}\n")    # messages
               
        # Line Analysis
        self.la_indir.set(self.nor_outdir.get())     # input folder
        if self.do_lsd.get():
            self.la_spec.set('*_lsd.fits')      # input (normalised) LSD profile OR pattern
        else:
            self.la_spec.set('*_ccf.fits')      # input (normalised) LSD profile OR pattern
        self.la_outdir.set(self.nor_outdir.get())     # input folder
        self.option_limits.set(self.option_la_values[0])  # choose the line limits definition
        self.show_limitup.set(u'3 \u03c3')     # show either manual limit or 3sigma/vsini
        self.show_limitlow.set(u'3 \u03c3')    # show either manual limit or 3sigma/vsini
        
        # Normalise LSD profiles 
        self.limitlow.set(-80)              # lower line limit (default self.rvmin + something)
        self.limitup.set(80)                # upper line limit (default self.rvmax - something)
        self.std.set(1)                     # create plot with mean profile and stdev
        self.norprfplot.set(1)              # plot the output
        self.norprf_output.set('_pfn.fits') # suffix of the output normalised LSD profiles (if missing, no output saved)
                
        # Fit LSD profiles
        self.fitgauss.set(1)                # Gaussian fit
        self.fitrot.set(1)                  # rotational fit
        self.fitlorentz.set(1)              # Lorentzian fit
        self.fitvoigt.set(1)                # Voigt fit
        self.rv0.set(0)                     # guess RV value
        self.width.set(10)                  # guess line width value
        self.ld.set(0.6)                    # linear limb darkening value
        self.fit_errs.set(1)                # use the data errors when fitting
        self.fitplot.set(1)                 # plot the output
        self.fit_output.set('_fit.fits')    # suffix of the fitting output (if missing, no output saved)
        
        # Compute moments
        self.usegauss.set(1)                # use Gaussian 3sigma to define the line limit
        self.userot.set(0)                  # use rotational width to define the line limit
        self.momplot.set(1)                 # plot the output
        self.mom_output.set('_mom.fits')    # suffix of the moments output (if missing, no output saved)
        
        # Compute bisector
        self.bisplot.set(1)                 # plot the output
        self.bis_output.set('_bis.fits')     # suffix of the bisector output (if missing, no output saved)
        
        # Compute Fourier transform
        self.fouplot.set(1)                 # plot the output
        self.fou_output.set('_fou.fits')     # suffix of the Fourier output (if missing, no output saved)
        
        # Threads
        self.thread_running.set(0)          # check if threads are running
        

    def define_plot(self):
        self.do_plot = False
        self.xvalues = None
        self.yvalues = None
        self.y_add = None
        self.y_res = None
        self.ymin = None
        self.ymax = None
        self.limits = (None,None)
        self.logscale = False
        self.close_plot = False
        self.plot_reset = False

    def create_widgets(self):

        # Define button style
        s = ttk.Style()
        s.configure('function.TButton', background='lightgreen')
        s.configure('abort.TButton', background='pink', font=(self.font_type, self.label_size, 'bold'))
        s.configure('browse.TButton', background='lightyellow')
        s.configure('function.TLabelframe.Label', font=(self.font_type, self.label_size, 'bold'))
        s.configure('subfunction.TLabelframe.Label', font=(self.font_type, self.text_size, 'bold'))

        # Create all empty frames
        
        # Normalisation frame (Upper left)
        #norrows = 6
        #norcols = 15
        norframe = ttk.LabelFrame(self.master, text='1. Normalisation', padding="3 3 3 3", style='function.TLabelframe')
        norframe.grid(row=0, rowspan=self.norrows, column=0, columnspan=self.norcols, sticky='nsew')        
        for n in range(self.norrows):
            norframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(self.norcols):
            norframe.columnconfigure(n, weight=1, minsize=1)
        
        # Profile frame (Lower left)
        #lsdrows = 11
        lsdcols = self.norcols        
        lsdframe = ttk.LabelFrame(self.master, text='2. Line Profile', padding="3 3 3 3", style='function.TLabelframe')
        lsdframe.grid(row=self.norrows, rowspan=self.lsdrows, column=0, columnspan=lsdcols, sticky='nsew')
        for n in range(self.lsdrows):
            lsdframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(lsdcols):
            lsdframe.columnconfigure(n, weight=1, minsize=1)

        # 1. Extract a single line (nested in Profile)
        extcols = self.norcols
        extframe = ttk.LabelFrame(lsdframe, text='2a. Single Line Extraction', padding="3 3 3 3", style='subfunction.TLabelframe')
        extframe.grid(row=self.lsdrows_input, rowspan=self.extrows, column=0, columnspan=extcols, sticky='nsew')
        for n in range(self.extrows):
            extframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(extcols):
            extframe.columnconfigure(n, weight=1, minsize=1)        
        
        # 2. Compute LSD/CCF Profile (nested in Profile)
        prfcols = self.norcols
        prfframe = ttk.LabelFrame(lsdframe, text='2b. Compute LSD or CCF', padding="3 3 3 3", style='subfunction.TLabelframe')
        prfframe.grid(row=self.lsdrows_input+self.extrows, rowspan=self.prfrows, column=0, columnspan=prfcols, sticky='nsew')
        for n in range(self.prfrows):
            prfframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(prfcols):
            prfframe.columnconfigure(n, weight=1, minsize=1) 
        
        # Plot frame (Upper centre)
        #self.pltrows = 7
        #self.pltcols = 11
        self.pltframe = ttk.LabelFrame(self.master, text='Plots', padding="3 3 3 3", style='function.TLabelframe')
        self.pltframe.grid(row=0, rowspan=self.pltrows, column=self.norcols, columnspan=self.pltcols, sticky='nsew')
        for n in range(self.pltrows):
            self.pltframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(self.pltcols):
            self.pltframe.columnconfigure(n, weight=1, minsize=1)
        
        # Message frame (Lower centre)
        msgcols = self.pltcols        
        msgframe = ttk.LabelFrame(self.master, text='Messages', padding="3 3 3 3", style='function.TLabelframe')
        msgframe.grid(row=self.pltrows, rowspan=self.msgrows, column=self.norcols, columnspan=msgcols, sticky='nsew')
        for n in range(self.msgrows):
            msgframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(msgcols):
            msgframe.columnconfigure(n, weight=1, minsize=1)
        
        # Line Analysis frame (Upper right)
        laframe = ttk.LabelFrame(self.master, text='3. Line Analysis', padding="3 3 3 3", style='function.TLabelframe')
        laframe.grid(row=0, rowspan=self.larows, column=self.norcols+self.pltcols, columnspan=self.lacols, sticky='nsew')
        for n in range(self.larows):
            laframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(self.lacols):
            laframe.columnconfigure(n, weight=1, minsize=1)
        
        # 1. Normalise Profile (nested in Line Analysis)
        norprfcols = self.lacols
        norprfframe = ttk.LabelFrame(laframe, text='3a. Profile Normalisation', padding="3 3 3 3", style='subfunction.TLabelframe')
        norprfframe.grid(row=self.larows_input, rowspan=self.norprfrows, column=0, columnspan=norprfcols, sticky='nsew')
        for n in range(self.norprfrows):
            norprfframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(norprfcols):
            norprfframe.columnconfigure(n, weight=1, minsize=1)        
        
        # 2. Fit Profile (nested in Line Analysis)
        fitcols = self.lacols
        fitframe = ttk.LabelFrame(laframe, text='3b. Profile Fitting', padding="3 3 3 3", style='subfunction.TLabelframe')
        fitframe.grid(row=self.larows_input+self.norprfrows, rowspan=self.fitrows, column=0, columnspan=fitcols, sticky='nsew')
        for n in range(self.fitrows):
            fitframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(fitcols):
            fitframe.columnconfigure(n, weight=1, minsize=1) 
        
        # 3. Moments (nested in Line Analysis)
        momcols = self.lacols
        momframe = ttk.LabelFrame(laframe, text='3c. Moments', padding="3 3 3 3", style='subfunction.TLabelframe')
        momframe.grid(row=self.larows_input+self.norprfrows+self.fitrows, rowspan=self.momrows, column=0, columnspan=momcols, sticky='nsew')
        self.option_mom_var = tk.StringVar(momframe)    # variable type for the OptionMenu
        for n in range(self.momrows):
            momframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(momcols):
            momframe.columnconfigure(n, weight=1, minsize=1) 
        
        # 4. Bisector (nested in Line Analysis)
        biscols = self.lacols
        bisframe = ttk.LabelFrame(laframe, text='3d. Bisector', padding="3 3 3 3", style='subfunction.TLabelframe')
        bisframe.grid(row=self.larows_input+self.norprfrows+self.fitrows+self.momrows, rowspan=self.bisrows, column=0, columnspan=biscols, sticky='nsew')
        self.option_bis_var = tk.StringVar(bisframe)    # variable type for the OptionMenu
        for n in range(self.bisrows):
            bisframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(biscols):
            bisframe.columnconfigure(n, weight=1, minsize=1) 
        
        # 5. Fourier transform (nested in Line Analysis)
        foucols = self.lacols        
        fouframe = ttk.LabelFrame(laframe, text='3e. Fourier Transform', padding="3 3 3 3", style='subfunction.TLabelframe')
        fouframe.grid(row=self.larows_input+self.norprfrows+self.fitrows+self.momrows+self.bisrows, rowspan=self.fourows, column=0, columnspan=foucols, sticky='nsew')
        self.option_fou_var = tk.StringVar(fouframe)    # variable type for the OptionMenu
        for n in range(self.fourows):
            fouframe.rowconfigure(n, weight=1, minsize=1)
        for n in range(foucols):
            fouframe.columnconfigure(n, weight=1, minsize=1) 


        # bottom frame (Lower right)
        bottom = ttk.Frame(self.master, padding="3 3 3 3")
        bottomcols = self.lacols
        bottom.grid(row=self.larows,rowspan=self.bottomrows, column=self.norcols+self.pltcols, columnspan=bottomcols, sticky='nsew')
        for n in range(self.bottomrows):
            bottom.rowconfigure(n, weight=1, minsize=1) #, uniform="a")
        for n in range(bottomcols):
            bottom.columnconfigure(n, weight=1, minsize=1) #, uniform="a")
            
            
        # Create all tab, field and buttons
        
        # 1. Normalisation frame
        # 1st row: select directory
        ttk.Label(norframe, text="Input folder: ").grid(row=0, column=0, columnspan=4, sticky='nsew')
        self.nor_indir_entry = ttk.Entry(norframe, textvariable=self.nor_indir, width=10)
        self.nor_indir_entry.grid(row=0, column=4, columnspan=10, sticky='nsew')
        self.nor_indir_entry.xview(tk.END)
        ttk.Button(norframe, text="Browse", command=self.nor_load_indir, style='browse.TButton').grid(row=0, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.nor_indir_entry, msg='Directory with input FITS/ASCII spectra', fg=self.t_tcolor, bg=self.t_background)

        # 2nd row: select file/pattern
        ttk.Label(norframe, text="File/Pattern: ").grid(row=1, column=0, columnspan=4, sticky='nsew')        
        self.nor_spec_entry = ttk.Entry(norframe, textvariable=self.nor_spec, width=10)
        self.nor_spec_entry.grid(row=1, column=4, columnspan=10, sticky='nsew')
        self.nor_spec_entry.xview(tk.END)
        ttk.Button(norframe, text="Browse", command=self.nor_load_file, style='browse.TButton').grid(row=1, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.nor_spec_entry, msg='Select pattern OR single FITS/ASCII spectrum', fg=self.t_tcolor, bg=self.t_background)

        # 3rd row: set data columns
        
        self.instruments = ttk.Combobox(norframe, state="readonly",\
                textvariable=self.option_instr,\
                values=self.option_values, width=12)
        self.instruments.bind("<<ComboboxSelected>>", self.change_cols)
        self.instruments.grid(row=2, column=0, columnspan=4, sticky='nsew')
        
        ToolTip(self.instruments, msg='Select an instrument to automatically define the columns/fields/hdus of the data.\nUNDEF may be used for:\n- ASCII files with wavelength and flux (2 columns)\n- FITS monodimensional file with the wavelength from CRVAL1, CDELT1, NAXIS1\n- all HARPS/HARPS-N/SOPHIE s1d and e2ds data (old DRS)\n  and GIANO-B s1d data', fg=self.t_tcolor, bg=self.t_background)
        
        ttk.Label(norframe, text="Wave:").grid(row=2, column=4, columnspan=2, sticky='nse')
        self.wavecol_entry = ttk.Entry(norframe, textvariable=self.wavecol, width=4)
        self.wavecol_entry.grid(row=2, column=6, sticky='nsw')
        ToolTip(self.wavecol_entry, msg='Column/Field with wavelength (1 = 1st column)\nRead the GIANO-B ms1d data with the options: Wave=2, Flux=3, S/N=4, Orders=1\nRead the ESPRESSO S1D data with the options: Wave=1, Flux=3, Errors=4\nRead the ESPRESSO S2D data with the options: Wave=4, Flux=1, Errors=2', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norframe, text="Flux:").grid(row=2, column=7, columnspan=2, sticky='nse')
        self.fluxcol_entry = ttk.Entry(norframe, textvariable=self.fluxcol, width=4)
        self.fluxcol_entry.grid(row=2, column=9, sticky='nsw')
        ToolTip(self.fluxcol_entry, msg='Column/Field with flux (1 = 1st column)\nRead the GIANO-B ms1d data with the options: Wave=2, Flux=3, S/N=4, Orders=1\nRead the ESPRESSO S1D data with the options: Wave=1, Flux=3, Errors=4\nRead the ESPRESSO S2D data with the options: Wave=4, Flux=1, Errors=2', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norframe, text="S/N:").grid(row=2, column=10, sticky='nse')
        self.snrcol_entry = ttk.Entry(norframe, textvariable=self.snrcol, width=4)
        self.snrcol_entry.grid(row=2, column=11, sticky='nsw')
        ToolTip(self.snrcol_entry, msg='Column/Field with S/N (1 = 1st column, 0 = no data)\nRead the GIANO-B ms1d data with the options: Wave=2, Flux=3, S/N=4, Orders=1\nRead the ESPRESSO S1D data with the options: Wave=1, Flux=3, Errors=4\nRead the ESPRESSO S2D data with the options: Wave=4, Flux=1, Errors=2', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norframe, text="Errs:").grid(row=2, column=12, sticky='nse')
        self.errcol_entry = ttk.Entry(norframe, textvariable=self.errcol, width=4)
        self.errcol_entry.grid(row=2, column=13, sticky='nsw')
        ToolTip(self.errcol_entry, msg='Column/Field with errors (1 = 1st column, 0 = no data)\nIf the S/N field is given), then the errors will not be used\nRead the GIANO-B ms1d data with the options: Wave=2, Flux=3, S/N=4, Orders=1\nRead the ESPRESSO S1D data with the options: Wave=1, Flux=3, Errors=4\nRead the ESPRESSO S2D data with the options: Wave=4, Flux=1, Errors=2', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norframe, text="Orders:").grid(row=2, column=14, columnspan=2, sticky='nse')
        self.echcol_entry = ttk.Entry(norframe, textvariable=self.echcol, width=4)
        self.echcol_entry.grid(row=2, column=16, sticky='nsw')
        ToolTip(self.echcol_entry, msg='Column/Field with echelle orders (1 = 1st column, 0 = no data)\nRead the GIANO-B ms1d data with the options: Wave=2, Flux=3, S/N=4, Orders=1\nRead the ESPRESSO S1D data with the options: Wave=1, Flux=3, Errors=4\nRead the ESPRESSO S2D data with the options: Wave=4, Flux=1, Errors=2', fg=self.t_tcolor, bg=self.t_background)


        # 4th row: wavelength units, fit parameters, plot
        self.units = ttk.Combobox(norframe, state="readonly",\
                textvariable=self.units_default,\
                values=self.units_values, width=10)
        self.units.bind("<<ComboboxSelected>>", self.change_units)
        self.units.grid(row=3, column=0, columnspan=4, sticky='nsew')
        ToolTip(self.units, msg='Select the wavelength unit', fg=self.t_tcolor, bg=self.t_background)

        self.spec_frame = ttk.Combobox(norframe, state="readonly",\
                textvariable=self.wave_frame,\
                values=self.wave_values, width=10)
        self.spec_frame.bind("<<ComboboxSelected>>", self.change_frame)        
        self.spec_frame.grid(row=3, column=4, columnspan=3, sticky='nsew')
        ToolTip(self.spec_frame, msg='Select the wavelength frame', fg=self.t_tcolor, bg=self.t_background)     
        
        ttk.Label(norframe, text="Degree:").grid(row=3, column=7, columnspan=2, sticky='nse')
        self.degree_entry = ttk.Entry(norframe, textvariable=self.degree, width=4)
        self.degree_entry.grid(row=3, column=9, columnspan=1, sticky='nsw')
        ToolTip(self.degree_entry, msg='Degree of polynomial for continuum fitting', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norframe, text="Subsets:").grid(row=3, column=10, columnspan=3, sticky='nse')
        self.n_ord_entry = ttk.Entry(norframe, textvariable=self.n_ord, width=5)
        self.n_ord_entry.grid(row=3, column=13, sticky='nsw')
        ToolTip(self.n_ord_entry, msg='ONLY without echelle orders given: normalise subsets independently', fg=self.t_tcolor, bg=self.t_background)
        #ttk.Label(norframe, text="Refine:", width=8).grid(row=3, column=10, sticky='nse')
        #self.refine_entry = ttk.Checkbutton(norframe, variable=self.refine, width=3)#, style="sp.Toolbutton", text="Refine")
        #self.refine_entry.grid(row=3, column=11, sticky='nsw')
        #ToolTip(self.refine_entry, msg='Refine normalisation', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norframe, text="Plot:").grid(row=3, column=14, columnspan=2, sticky='nse')
        self.norplot_entry = ttk.Checkbutton(norframe, variable=self.norplot, width=3)
        self.norplot_entry.grid(row=3, column=16, sticky='nsw')
        ToolTip(self.norplot_entry, msg='Plot normalised spectra', fg=self.t_tcolor, bg=self.t_background)
        
        # 5th row: output directory
        ttk.Label(norframe, text="Output folder: ").grid(row=4, column=0, columnspan=4, sticky='nsew')
        self.nor_outdir_entry = ttk.Entry(norframe, textvariable=self.nor_outdir, width=10)
        self.nor_outdir_entry.grid(row=4, column=4, columnspan=10, sticky='nsew')
        self.nor_outdir_entry.xview(tk.END)
        ttk.Button(norframe, text="Browse", command=self.nor_load_outdir, style='browse.TButton').grid(row=4, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.nor_outdir_entry, msg='Output directory', fg=self.t_tcolor, bg=self.t_background)
        
        # 6th row: output suffix and run normalisation
        ttk.Label(norframe, text="Suffix: ").grid(row=5, column=0, columnspan=4, sticky='nsew')
        self.nor_output_entry = ttk.Entry(norframe, textvariable=self.nor_output, width=10)
        self.nor_output_entry.grid(row=5, column=4, columnspan=6, sticky='nsew')
        self.nor_output_entry.xview(tk.END)
        ToolTip(self.nor_output_entry, msg='Suffix for normalised FITS spectra', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(norframe, text="Normalise spectra", style='function.TButton', command=self.normalise).grid(row=5, column=10, columnspan=7, sticky='nsew')

        for child in norframe.winfo_children(): child.grid_configure(padx=5, pady=5)
        
        # 2. Profile frame
        # Input/output folder and data
        # 1st row: select directory
        ttk.Label(lsdframe, text="Input folder: ").grid(row=0, column=0, columnspan=4, sticky='nsew')
        self.lsd_indir_entry = ttk.Entry(lsdframe, textvariable=self.lsd_indir, width=10)
        self.lsd_indir_entry.grid(row=0, column=4, columnspan=10, sticky='nsew')
        self.lsd_indir_entry.xview(tk.END)
        ttk.Button(lsdframe, text="Browse", command=self.lsd_load_indir, style='browse.TButton').grid(row=0, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.lsd_indir_entry, msg='Directory with input normalised spectra', fg=self.t_tcolor, bg=self.t_background)
        
        # 2nd row: select file/pattern
        ttk.Label(lsdframe, text="File/Pattern: ").grid(row=1, column=0, columnspan=4, sticky='nsew')
        self.lsd_spec_entry = ttk.Entry(lsdframe, textvariable=self.lsd_spec, width=10)
        self.lsd_spec_entry.grid(row=1, column=4, columnspan=10, sticky='nsew')
        self.lsd_spec_entry.xview(tk.END)
        ttk.Button(lsdframe, text="Browse", command=self.lsd_load_file, style='browse.TButton').grid(row=1, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.lsd_spec_entry, msg='Select pattern OR single normalised spectrum\nThe spectrum MUST be the output of the normalisation procedure:\nwavelength in vacuum and in angstroms\nFields: wave=1, flux=2, Nflux=3, S/N=4', fg=self.t_tcolor, bg=self.t_background)  
            
        
        # 3rd row: output directory
        ttk.Label(lsdframe, text="Output folder: ").grid(row=2, column=0, columnspan=4, sticky='nsew')
        self.lsd_outdir_entry = ttk.Entry(lsdframe, textvariable=self.lsd_outdir, width=10)
        self.lsd_outdir_entry.grid(row=2, column=4, columnspan=10, sticky='nsew')
        self.lsd_outdir_entry.xview(tk.END)
        ttk.Button(lsdframe, text="Browse", command=self.lsd_load_outdir, style='browse.TButton').grid(row=2, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.lsd_outdir_entry, msg='Output directory', fg=self.t_tcolor, bg=self.t_background)
        
        # 2a. Extract single line
        # 1st row: define wavelength, RV range, plot      
        ttk.Label(extframe, text="Wavelength:").grid(row=0, column=0, columnspan=4, sticky='nsw')
        self.extwave_entry = ttk.Entry(extframe, textvariable=self.ext_wave, width=10)
        self.extwave_entry.grid(row=0, column=4, columnspan=2, sticky='nsw')      
        ToolTip(self.extwave_entry, msg='Wavelength of the line to extract, in Angstroms. It will define RV=0', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(extframe, text="RV min: ").grid(row=0, column=6, columnspan=2, sticky='nse')
        self.rvmin_entry = ttk.Entry(extframe, textvariable=self.rvmin, width=5)
        self.rvmin_entry.grid(row=0, column=8, columnspan=1, sticky='nsw')      
        ToolTip(self.rvmin_entry, msg='Minimum RV of the profile', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(extframe, text="RV max: ").grid(row=0, column=9, columnspan=2, sticky='nse')
        self.rvmax_entry = ttk.Entry(extframe, textvariable=self.rvmax, width=5)
        self.rvmax_entry.grid(row=0, column=11, columnspan=1, sticky='nsw')       
        ToolTip(self.rvmax_entry, msg='Maximum RV of the profile', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(extframe, text="Step: ").grid(row=0, column=12, columnspan=2, sticky='nse')
        self.rvstep_entry = ttk.Entry(extframe, textvariable=self.rvstep, width=3)
        self.rvstep_entry.grid(row=0, column=14, columnspan=1, sticky='nsw')
        ToolTip(self.rvstep_entry, msg='RV step of the profile', fg=self.t_tcolor, bg=self.t_background)      
        ttk.Label(extframe, text="Plot: ").grid(row=0, column=15, columnspan=1, sticky='nse')
        self.extplot_entry = ttk.Checkbutton(extframe, variable=self.extplot, width=3)
        self.extplot_entry.grid(row=0, column=16, sticky='nsw')
        ToolTip(self.extplot_entry, msg='Plot profiles', fg=self.t_tcolor, bg=self.t_background)
        
        # 2nd row: output suffix and run extraction   
        ttk.Label(extframe, text="Suffix: ").grid(row=1, column=0, columnspan=4, sticky='nsew')
        self.ext_output_entry = ttk.Entry(extframe, textvariable=self.ext_output, width=10)
        self.ext_output_entry.grid(row=1, column=4, columnspan=6, sticky='nsew')
        self.ext_output_entry.xview(tk.END)
        ToolTip(self.ext_output_entry, msg='Suffix for single line output profiles', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(extframe, text="Extract line", style='function.TButton', command=self.thread_extract_line).grid(row=1, column=10, columnspan=7, sticky='nsew')

        for child in extframe.winfo_children(): child.grid_configure(padx=5, pady=5)        

        # 2b. Compute LSD or CCF profile
        # 1st row: select mask
        ttk.Label(prfframe, text="Mask: ").grid(row=0, column=0, columnspan=3, sticky='nsew')
        self.mask_entry = ttk.Entry(prfframe, textvariable=self.mask, width=10)
        self.mask_entry.grid(row=0, column=3, columnspan=5, sticky='nsew')
        self.mask_entry.xview(tk.END)
        ttk.Button(prfframe, text="Browse", command=self.mask_load_file, style='browse.TButton').grid(row=0, column=8, columnspan=3, sticky='nsew')
        ToolTip(self.mask_entry, msg='Select mask (VALD file or 2-column ASCII file)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="Invert: ").grid(row=0, column=11, columnspan=2, sticky='nse')
        self.invert_entry = ttk.Checkbutton(prfframe, variable=self.mask_invert, width=3)
        self.invert_entry.grid(row=0, column=13, sticky='nsw')
        ToolTip(self.invert_entry, msg='Convert fluxes to depths (absorption mask/spectrum/model ONLY)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="Spectrum: ").grid(row=0, column=14, columnspan=2, sticky='nse')
        self.maskspec_entry = ttk.Checkbutton(prfframe, variable=self.mask_spectrum, width=3)
        self.maskspec_entry.grid(row=0, column=16, sticky='nsw')
        ToolTip(self.maskspec_entry, msg='Use a spectrum (with continuum) as mask', fg=self.t_tcolor, bg=self.t_background)
        # 2nd row: units, depths, wavelength range
        self.mask_units = ttk.Combobox(prfframe, state="readonly",\
                textvariable=self.mask_units_default,\
                values=self.mask_units_values, width=10)
        self.mask_units.bind("<<ComboboxSelected>>", self.mask_change_units)  
        self.mask_units.grid(row=1, column=0, columnspan=3, sticky='nsew')
        ToolTip(self.mask_units, msg='Select the wavelength unit', fg=self.t_tcolor, bg=self.t_background)
        
        self.mask_frame = ttk.Combobox(prfframe, state="readonly",\
                textvariable=self.mask_wave_frame,\
                values=self.mask_wave_values, width=10)
        self.mask_frame.bind("<<ComboboxSelected>>", self.mask_change_frame) 
        self.mask_frame.grid(row=1, column=3, columnspan=2, sticky='nsew')
        ToolTip(self.mask_frame, msg='Select the wavelength frame', fg=self.t_tcolor, bg=self.t_background)   
        
        ttk.Label(prfframe, text="D min: ").grid(row=1, column=5, columnspan=2, sticky='nse')
        self.mask_dlow_entry = ttk.Entry(prfframe, textvariable=self.mask_dlow, width=5)
        self.mask_dlow_entry.grid(row=1, column=7, sticky='nsw')        
        ToolTip(self.mask_dlow_entry, msg='Minimum line depth of mask lines to be used', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="D max: ").grid(row=1, column=8, columnspan=2, sticky='nse')
        self.mask_dup_entry = ttk.Entry(prfframe, textvariable=self.mask_dup, width=5)
        self.mask_dup_entry.grid(row=1, column=10, sticky='nsw')        
        ToolTip(self.mask_dup_entry, msg='Maximum line depth of mask lines to be used', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="W min: ").grid(row=1, column=11, columnspan=2, sticky='nse')
        self.mask_wmin_entry = ttk.Entry(prfframe, textvariable=self.mask_wmin, width=5)
        self.mask_wmin_entry.grid(row=1, column=13, columnspan=1, sticky='nsw')     
        ToolTip(self.mask_wmin_entry, msg='Minimum wavelength of mask lines to be used', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="W max: ").grid(row=1, column=14, columnspan=2, sticky='nse')
        self.mask_wmax_entry = ttk.Entry(prfframe, textvariable=self.mask_wmax, width=5)
        self.mask_wmax_entry.grid(row=1, column=16, columnspan=1, sticky='nsw')      
        ToolTip(self.mask_wmax_entry, msg='Maximum wavelength of mask lines to be used', fg=self.t_tcolor, bg=self.t_background)
        
        # 3rd row: select elements
        ttk.Label(prfframe, text="VALD els: ").grid(row=2, column=0, columnspan=3, sticky='nsew')
        self.mask_els_entry = ttk.Entry(prfframe, textvariable=self.mask_els, width=15)
        self.mask_els_entry.grid(row=2, column=3, columnspan=5, sticky='nsew')
        ToolTip(self.mask_els_entry, msg='ONLY for VALD mask: use only selected elements, e.g. "Fe 1,Fe 2,H 1"', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="VALD NO els: ").grid(row=2, column=8, columnspan=4, sticky='nse')
        self.mask_noels_entry = ttk.Entry(prfframe, textvariable=self.mask_noels, width=15)
        self.mask_noels_entry.grid(row=2, column=12, columnspan=5, sticky='nsew')
        ToolTip(self.mask_noels_entry, msg='ONLY for VALD mask: exclude selected elements, e.g. "Fe 1,Fe 2,H 1"', fg=self.t_tcolor, bg=self.t_background)
        
        # 4th row: balmer, tellurics, remove cosmic rays, clean, S/N as weights
        ttk.Label(prfframe, text="Balmer: ", width=8).grid(row=3, column=0, columnspan=2, sticky='nse')
        self.mask_balmer_entry = ttk.Checkbutton(prfframe, variable=self.mask_balmer, width=3)
        self.mask_balmer_entry.grid(row=3, column=2, sticky='nsw')
        ToolTip(self.mask_balmer_entry, msg='Exclude Balmer regions', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="Tellurics: ").grid(row=3, column=3, columnspan=3, sticky='nse')
        self.mask_tell_entry = ttk.Checkbutton(prfframe, variable=self.mask_tell, width=3)
        self.mask_tell_entry.grid(row=3, column=6, sticky='nsw')
        ToolTip(self.mask_tell_entry, msg='Exclude telluric regions', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="Cosmics: ").grid(row=3, column=7, columnspan=2, sticky='nse')
        self.cosmic_entry = ttk.Checkbutton(prfframe, variable=self.cosmic, width=3)
        self.cosmic_entry.grid(row=3, column=9, sticky='nsw')
        ToolTip(self.cosmic_entry, msg='Remove cosmic rays', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="Clean: ").grid(row=3, column=10, columnspan=2, sticky='nse')
        self.clean_entry = ttk.Checkbutton(prfframe, variable=self.clean, width=3)
        self.clean_entry.grid(row=3, column=12, sticky='nsw')
        ToolTip(self.clean_entry, msg='Clean spectrum (smoothing)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="S/N weighted: ").grid(row=3, column=13, columnspan=3, sticky='nse')
        self.ccfweight_entry = ttk.Checkbutton(prfframe, variable=self.ccfweight, width=3)
        self.ccfweight_entry.grid(row=3, column=16, sticky='nsw')
        ToolTip(self.ccfweight_entry, msg='Use the normalised S/N values as weigths when computing the CCF', fg=self.t_tcolor, bg=self.t_background) 
        
        # 5th row: rv range, rv step, LSD, CCF, plot
        
        ttk.Label(prfframe, text="RV min: ").grid(row=4, column=0, columnspan=2, sticky='nse')
        self.rvmin_entry = ttk.Entry(prfframe, textvariable=self.rvmin, width=5)
        self.rvmin_entry.grid(row=4, column=2, columnspan=2, sticky='nsw')      
        ToolTip(self.rvmin_entry, msg='Minimum RV of the profile', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="RV max: ", width=8).grid(row=4, column=4, columnspan=2, sticky='nse')
        self.rvmax_entry = ttk.Entry(prfframe, textvariable=self.rvmax, width=5)
        self.rvmax_entry.grid(row=4, column=6, columnspan=2, sticky='nsw')       
        ToolTip(self.rvmax_entry, msg='Maximum RV of the profile', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="Step: ").grid(row=4, column=8, columnspan=1, sticky='nse')
        self.rvstep_entry = ttk.Entry(prfframe, textvariable=self.rvstep, width=3)
        self.rvstep_entry.grid(row=4, column=9, columnspan=1, sticky='nsw')
        ToolTip(self.rvstep_entry, msg='RV step of the profile', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="LSD: ").grid(row=4, column=11, columnspan=1, sticky='nse')
        self.lsdlsd_entry = ttk.Checkbutton(prfframe, variable=self.do_lsd, command=self.sfx_lsd, width=3)
        self.lsdlsd_entry.grid(row=4, column=12, sticky='nsw')
        ToolTip(self.lsdlsd_entry, msg='Compute LSD profiles', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="CCF: ").grid(row=4, column=13, columnspan=1, sticky='nse')
        self.lsdccf_entry = ttk.Checkbutton(prfframe, variable=self.do_ccf, command=self.sfx_ccf, width=3)
        self.lsdccf_entry.grid(row=4, column=14, sticky='nsw')
        ToolTip(self.lsdccf_entry, msg='Compute CCF profiles', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(prfframe, text="Plot: ").grid(row=4, column=15, columnspan=1, sticky='nse')
        self.lsdplot_entry = ttk.Checkbutton(prfframe, variable=self.lsdplot, width=3)
        self.lsdplot_entry.grid(row=4, column=16, sticky='nsw')
        ToolTip(self.lsdplot_entry, msg='Plot profiles', fg=self.t_tcolor, bg=self.t_background)   
          
        # 6th row: output suffix and run LSD/CCF
        ttk.Label(prfframe, text="Suffix: ").grid(row=5, column=0, columnspan=3, sticky='nsew')
        self.lsd_output_entry = ttk.Entry(prfframe, textvariable=self.lsd_output, width=10)
        self.lsd_output_entry.grid(row=5, column=3, columnspan=6, sticky='nsew')
        self.lsd_output_entry.xview(tk.END)
        ToolTip(self.lsd_output_entry, msg='Suffix for FITS output profiles', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(prfframe, text="Compute profile", style='function.TButton', command=self.thread_do_profile).grid(row=5, column=9, columnspan=8, sticky='nsew')

        for child in prfframe.winfo_children(): child.grid_configure(padx=5, pady=5)
        for child in lsdframe.winfo_children(): child.grid_configure(padx=5, pady=5)
        
        
        # PLOT FRAME        

        self.thread_plot()
        #self.fig, self.ax = plt.subplots(figsize=(4,3), dpi=100)
        #self.plot_logo()

        
        # MESSAGE FRAME
        
        self.messages_entry = tk.Text(msgframe, height=self.msgrows,width=40)
        self.messages_entry.grid(row=0, rowspan=self.msgrows-1, column=0, columnspan=msgcols-1, sticky='nsew')
        
        # create a Scrollbar and associate it with txt
        self.scrollbar = ttk.Scrollbar(msgframe, command=self.messages_entry.yview)
        self.scrollbar.grid(row=0, rowspan=self.msgrows-1, column=msgcols-1, sticky='nsew')
        self.messages_entry['yscrollcommand'] = self.scrollbar.set
        
        # Custom font
        cf = tk.font.Font(family=self.font_type, size=self.text_size)
        self.messages_entry.configure(font=cf)
        self.messages_entry.insert(tk.END, self.messages.get())
        
        ttk.Button(msgframe, text="Save Log", style='function.TButton', command=self.save_log).grid(row=self.msgrows-1, column=int(msgcols/2)-1, columnspan=3, sticky='nsew')
        
        for child in msgframe.winfo_children(): child.grid_configure(padx=5, pady=5)
        
        
        # 3. Line Analysis frame
        # 3. Line Analysis
        # 3 - 1st row: select directory
        ttk.Label(laframe, text="Input folder: ").grid(row=0, column=0, columnspan=4, sticky='nsew')
        self.la_indir_entry = ttk.Entry(laframe, textvariable=self.la_indir, width=10)
        self.la_indir_entry.grid(row=0, column=4, columnspan=10, sticky='nsew')
        self.la_indir_entry.xview(tk.END)
        ttk.Button(laframe, text="Browse", command=self.la_load_indir, style='browse.TButton').grid(row=0, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.la_indir_entry, msg='Directory with input (normalised) profiles', fg=self.t_tcolor, bg=self.t_background)
        # 3 - 2nd row: select file/pattern
        ttk.Label(laframe, text="File/Pattern: ").grid(row=1, column=0, columnspan=4, sticky='nsew')
        self.la_spec_entry = ttk.Entry(laframe, textvariable=self.la_spec, width=10)
        self.la_spec_entry.grid(row=1, column=4, columnspan=10, sticky='nsew')
        self.la_spec_entry.xview(tk.END)
        ttk.Button(laframe, text="Browse", command=self.la_load_file, style='browse.TButton').grid(row=1, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.la_spec_entry, msg='Select pattern OR single (normalised) profiles', fg=self.t_tcolor, bg=self.t_background)
        # 3 - 3rd row: select output directory
        ttk.Label(laframe, text="Output folder: ").grid(row=2, column=0, columnspan=4, sticky='nsew')
        self.la_outdir_entry = ttk.Entry(laframe, textvariable=self.la_outdir, width=10)
        self.la_outdir_entry.grid(row=2, column=4, columnspan=10, sticky='nsew')
        self.la_outdir_entry.xview(tk.END)
        ttk.Button(laframe, text="Browse", command=self.la_load_outdir, style='browse.TButton').grid(row=2, column=14, columnspan=3, sticky='nsew')
        ToolTip(self.la_outdir_entry, msg='Output directory', fg=self.t_tcolor, bg=self.t_background)
        # 3a. Normalise profiles  
        # 3a - 1st row: limits, st. dev, plot
        ttk.Label(norprfframe, text="RV min: ").grid(row=0, column=0, columnspan=3, sticky='nse')
        self.limitlow_entry = ttk.Entry(norprfframe, textvariable=self.limitlow, width=3)
        self.limitlow_entry.grid(row=0, column=3, columnspan=2, sticky='nsw')
        ToolTip(self.limitlow_entry, msg='Left border of the profile (for continuum definition)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norprfframe, text="RV max: ").grid(row=0, column=5, columnspan=3, sticky='nse')
        self.limitup_entry = ttk.Entry(norprfframe, textvariable=self.limitup, width=3)
        self.limitup_entry.grid(row=0, column=8, columnspan=2, sticky='nsw')       
        ToolTip(self.limitup_entry, msg='Right border of the profile (for continuum definition)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norprfframe, text="St.Dev: ").grid(row=0, column=10, columnspan=3, sticky='nse')
        self.std_entry = ttk.Checkbutton(norprfframe, variable=self.std, width=3)
        self.std_entry.grid(row=0, column=13, sticky='nsw')
        ToolTip(self.std_entry, msg='Compute average and st. dev.', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(norprfframe, text="Plot: ").grid(row=0, column=14, columnspan=2, sticky='nse')
        self.norprfplot_entry = ttk.Checkbutton(norprfframe, variable=self.norprfplot, width=3)
        self.norprfplot_entry.grid(row=0, column=16, sticky='nsw')
        ToolTip(self.norprfplot_entry, msg='Plot normalised profiles', fg=self.t_tcolor, bg=self.t_background)  
        # 3a - 2nd row: suffix and run profile normalisation
        ttk.Label(norprfframe, text="Suffix: ").grid(row=1, column=0, columnspan=3, sticky='nsew')
        self.norprf_output_entry = ttk.Entry(norprfframe, textvariable=self.norprf_output, width=10)
        self.norprf_output_entry.grid(row=1, column=3, columnspan=7, sticky='nsew')
        ToolTip(self.norprf_output_entry, msg='Suffix for normalised LSD profiles', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(norprfframe, text="Normalise profile", style='function.TButton', command=self.thread_norm_profile).grid(row=1, column=10, columnspan=7, sticky='nsew')
        for child in norprfframe.winfo_children(): child.grid_configure(padx=5, pady=5)       
        # 3b. Fit profiles  
        # 3b - 1st row: fit function
        ttk.Label(fitframe, text="Gaussian: ", width=10).grid(row=0, column=0, columnspan=3, sticky='nse')
        self.gauss_entry = ttk.Checkbutton(fitframe, variable=self.fitgauss, width=3)
        self.gauss_entry.grid(row=0, column=3, sticky='nsw')
        ToolTip(self.gauss_entry, msg='Fit profile with a Gaussian function', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(fitframe, text="Lorentzian: ").grid(row=0, column=4, columnspan=4, sticky='nse')
        self.lor_entry = ttk.Checkbutton(fitframe, variable=self.fitlorentz, width=3)
        self.lor_entry.grid(row=0, column=8, sticky='nsw')       
        ToolTip(self.lor_entry, msg='Fit profile with a Lorentzian function', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(fitframe, text="Voigt: ").grid(row=0, column=9, columnspan=2, sticky='nse')
        self.voigt_entry = ttk.Checkbutton(fitframe, variable=self.fitvoigt, width=3)
        self.voigt_entry.grid(row=0, column=11, sticky='nsw')       
        ToolTip(self.voigt_entry, msg='Fit profile with a Voigt function', fg=self.t_tcolor, bg=self.t_background)  
        ttk.Label(fitframe, text="Rotational: ").grid(row=0, column=12, columnspan=4, sticky='nse')
        self.rot_entry = ttk.Checkbutton(fitframe, variable=self.fitrot, width=3)
        self.rot_entry.grid(row=0, column=16, sticky='nsw')       
        ToolTip(self.rot_entry, msg='Fit profile with a rotational function', fg=self.t_tcolor, bg=self.t_background)
        # 3b - 2nd row: plot parameters and plot
        ttk.Label(fitframe, text="RV guess: ", width=10).grid(row=1, column=0, columnspan=3, sticky='nse')
        self.rv0_entry = ttk.Entry(fitframe, textvariable=self.rv0, width=5)
        self.rv0_entry.grid(row=1, column=3, columnspan=1, sticky='nsw')
        ToolTip(self.rv0_entry, msg='RV guess value', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(fitframe, text="W guess: ", width=10).grid(row=1, column=4, columnspan=3, sticky='nse')
        self.width_entry = ttk.Entry(fitframe, textvariable=self.width, width=5)
        self.width_entry.grid(row=1, column=7, columnspan=1, sticky='nsw')
        ToolTip(self.width_entry, msg='Line width guess value', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(fitframe, text="LD: ", width=4).grid(row=1, column=8, columnspan=1, sticky='nse')
        self.ld_entry = ttk.Entry(fitframe, textvariable=self.ld, width=5)
        self.ld_entry.grid(row=1, column=9, sticky='nsw')        
        ToolTip(self.ld_entry, msg='Linear limb darkening  value', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(fitframe, text="Use errors: ", width=12).grid(row=1, column=10, columnspan=3, sticky='nse')
        self.fit_errs_entry = ttk.Checkbutton(fitframe, variable=self.fit_errs, width=3)
        self.fit_errs_entry.grid(row=1, column=13, sticky='nsw')        
        ToolTip(self.fit_errs_entry, msg='Use the errors on the data when fitting', fg=self.t_tcolor, bg=self.t_background)    
        ttk.Label(fitframe, text="Plot: ", width=6).grid(row=1, column=14, columnspan=2, sticky='nse')
        self.fitplot_entry = ttk.Checkbutton(fitframe, variable=self.fitplot, width=3)
        self.fitplot_entry.grid(row=1, column=16, sticky='nsw')
        ToolTip(self.fitplot_entry, msg='Plot fitting results', fg=self.t_tcolor, bg=self.t_background)
        
        # 3b - 3rd row: suffix and run profile fitting
        ttk.Label(fitframe, text="Suffix: ").grid(row=2, column=0, columnspan=3, sticky='nsew')
        self.fit_output_entry = ttk.Entry(fitframe, textvariable=self.fit_output, width=10)
        self.fit_output_entry.grid(row=2, column=3, columnspan=7, sticky='nsew')
        ToolTip(self.fit_output_entry, msg='Suffix for the fitting results', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(fitframe, text="Fit profile", style='function.TButton', command=self.thread_fit_profile).grid(row=2, column=10, columnspan=7, sticky='nsew')
        
        for child in fitframe.winfo_children(): child.grid_configure(padx=5, pady=5)
        
        # 3c. Compute Moments
        # 3c - 1st row: define limits and plot
        #self.mom_menu = ttk.OptionMenu(momframe, self.option_mom_var, self.option_la_default, *self.option_la_values, command=self.option_la_change)
        
        self.mom_menu = ttk.Combobox(momframe, state="readonly",\
                textvariable=self.option_limits,\
                values=self.option_la_values)
        self.mom_menu.bind("<<ComboboxSelected>>", self.option_la_change)
        self.mom_menu.grid(row=0, column=0, columnspan=4, sticky='nsew')
        ToolTip(self.mom_menu, msg='Select the line limits (from previous fit or manual)', fg=self.t_tcolor, bg=self.t_background)   
        
        ttk.Label(momframe, text="RV min: ").grid(row=0, column=4, columnspan=3, sticky='nse')
        self.limitlow_entry = ttk.Entry(momframe, textvariable=self.show_limitlow, width=6)
        self.limitlow_entry.grid(row=0, column=7, columnspan=2, sticky='nsw')        
        ToolTip(self.limitlow_entry, msg='Left border of the profile (for line definition)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(momframe, text="RV max: ").grid(row=0, column=9, columnspan=3, sticky='nse')
        self.limitup_entry = ttk.Entry(momframe, textvariable=self.show_limitup, width=6)
        self.limitup_entry.grid(row=0, column=12, columnspan=2, sticky='nsw')        
        ToolTip(self.limitup_entry, msg='Right border of the profile (for line definition)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(momframe, text="Plot: ").grid(row=0, column=14, columnspan=2, sticky='nse')
        self.momplot_entry = ttk.Checkbutton(momframe, variable=self.momplot, width=3)
        self.momplot_entry.grid(row=0, column=16, sticky='nsw')
        ToolTip(self.momplot_entry, msg='Plot line limits', fg=self.t_tcolor, bg=self.t_background)          
        # 3c - 2nd row: output suffix and run routine
        ttk.Label(momframe, text="Suffix: ").grid(row=1, column=0, columnspan=3, sticky='nsew')
        self.mom_output_entry = ttk.Entry(momframe, textvariable=self.mom_output, width=10)
        self.mom_output_entry.grid(row=1, column=3, columnspan=7, sticky='nsew')
        ToolTip(self.mom_output_entry, msg='Suffix for the moments results', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(momframe, text="Moments", style='function.TButton', command=self.thread_moments).grid(row=1, column=10, columnspan=7, sticky='nsew')
        
        for child in momframe.winfo_children(): child.grid_configure(padx=5, pady=5)
        
        # 3d. Compute Bisector
        # 3d - 1st row: define limits, plot
        #self.bis_menu = ttk.OptionMenu(bisframe, self.option_bis_var, self.option_la_default, *self.option_la_values, command=self.option_la_change)
        
        self.bis_menu = ttk.Combobox(bisframe, state="readonly",\
                textvariable=self.option_limits,\
                values=self.option_la_values)
        self.bis_menu.bind("<<ComboboxSelected>>", self.option_la_change)
        self.bis_menu.grid(row=0, column=0, columnspan=4, sticky='nsew')
        ToolTip(self.bis_menu, msg='Select the line limits (from previous fit or manual)', fg=self.t_tcolor, bg=self.t_background)  
        
        ttk.Label(bisframe, text="RV min: ").grid(row=0, column=4, columnspan=3, sticky='nse')
        self.limitlow_entry = ttk.Entry(bisframe, textvariable=self.show_limitlow, width=6)
        self.limitlow_entry.grid(row=0, column=7, columnspan=2, sticky='nsw')        
        ToolTip(self.limitlow_entry, msg='Left border of the profile (for line definition)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(bisframe, text="RV max: ").grid(row=0, column=9, columnspan=3, sticky='nse')
        self.limitup_entry = ttk.Entry(bisframe, textvariable=self.show_limitup, width=6)
        self.limitup_entry.grid(row=0, column=12, columnspan=2, sticky='nsw')        
        ToolTip(self.limitup_entry, msg='Right border of the profile (for line definition)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(bisframe, text="Plot: ").grid(row=0, column=14, columnspan=2, sticky='nse')
        self.bisplot_entry = ttk.Checkbutton(bisframe, variable=self.bisplot, width=3)
        self.bisplot_entry.grid(row=0, column=16, sticky='nsw')
        ToolTip(self.bisplot_entry, msg='Plot bisectors', fg=self.t_tcolor, bg=self.t_background)
        # 3d - 2nd row: output suffix and run routine
        ttk.Label(bisframe, text="Suffix: ").grid(row=1, column=0, columnspan=3, sticky='nsew')
        self.bis_output_entry = ttk.Entry(bisframe, textvariable=self.bis_output, width=10)
        self.bis_output_entry.grid(row=1, column=3, columnspan=7, sticky='nsew')
        ToolTip(self.bis_output_entry, msg='Suffix for the bisector results', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(bisframe, text="Bisector", style='function.TButton', command=self.thread_bisector).grid(row=1, column=10, columnspan=7, sticky='nsew')
        
        for child in bisframe.winfo_children(): child.grid_configure(padx=5, pady=5)
        
        # 3e. Compute Fourier Transform
        # 3e - 1st row: define limits
        #self.fou_menu = ttk.OptionMenu(fouframe, self.option_fou_var, self.option_la_default, *self.option_la_values, command=self.option_la_change)
        
        self.fou_menu = ttk.Combobox(fouframe, state="readonly",\
                textvariable=self.option_limits,\
                values=self.option_la_values)
        self.fou_menu.bind("<<ComboboxSelected>>", self.option_la_change)
        self.fou_menu.grid(row=0, column=0, columnspan=4, sticky='nsew')
        ToolTip(self.fou_menu, msg='Select the line limits (from previous fit or manual)', fg=self.t_tcolor, bg=self.t_background) 
        
        ttk.Label(fouframe, text="RV min: ").grid(row=0, column=4, columnspan=3, sticky='nse')
        self.limitlow_entry = ttk.Entry(fouframe, textvariable=self.show_limitlow, width=6)
        self.limitlow_entry.grid(row=0, column=7, columnspan=2, sticky='nsw')        
        ToolTip(self.limitlow_entry, msg='Left border of the profile (for line definition)', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(fouframe, text="RV max: ").grid(row=0, column=9, columnspan=3, sticky='nse')
        self.limitup_entry = ttk.Entry(fouframe, textvariable=self.show_limitup, width=6)
        self.limitup_entry.grid(row=0, column=12, columnspan=2, sticky='nsw')        
        ToolTip(self.limitup_entry, msg='Right border of the profile (for line definition)', fg=self.t_tcolor, bg=self.t_background)       
        ttk.Label(fouframe, text="Plot: ").grid(row=0, column=14, columnspan=2, sticky='nse')
        self.fouplot_entry = ttk.Checkbutton(fouframe, variable=self.fouplot, width=3)
        self.fouplot_entry.grid(row=0, column=16, sticky='nsw')
        ToolTip(self.fouplot_entry, msg='Plot Fourier Transform', fg=self.t_tcolor, bg=self.t_background)
        # 3e - 2nd row: limb darkening, output suffix and run routine
        ttk.Label(fouframe, text="LD: ").grid(row=1, column=0, sticky='nse')
        self.ld_entry = ttk.Entry(fouframe, textvariable=self.ld, width=3)
        self.ld_entry.grid(row=1, column=1, sticky='nsw')
        ToolTip(self.ld_entry, msg='Linear limb darkening value', fg=self.t_tcolor, bg=self.t_background)
        ttk.Label(fouframe, text="Suffix: ").grid(row=1, column=2, columnspan=3, sticky='nsew')
        self.fou_output_entry = ttk.Entry(fouframe, textvariable=self.fou_output, width=10)
        self.fou_output_entry.grid(row=1, column=5, columnspan=5, sticky='nsew')
        ToolTip(self.fou_output_entry, msg='Suffix for the Fourier Transform results', fg=self.t_tcolor, bg=self.t_background)
        ttk.Button(fouframe, text="Fourier Transform", style='function.TButton', command=self.fourier).grid(row=1, column=10, columnspan=7, sticky='nsew')

        
        for child in fouframe.winfo_children(): child.grid_configure(padx=5, pady=5)
 
 
        # bottom frame

        ttk.Button(bottom, text="Reset all fields", style='function.TButton', command=self.reset).grid(row=1, column=2, columnspan=4, sticky='nsew')
        
        ttk.Button(bottom, text="ABORT", style='abort.TButton', command=self.abort).grid(row=1, column=7, columnspan=4, sticky='nsew')

        ttk.Button(bottom, text="Quit", style='function.TButton', command=self.quit_shiva).grid(row=1, column=12, columnspan=4, sticky='nsew')
        
        ttk.Label(bottom, text=f"{self.version} -- {self.pa_version}").grid(row=2, column=2, columnspan=self.lacols-4, sticky='nsew')
        ttk.Label(bottom, text=f"GNU GPL v3 -- Monica Rainer -- INAF-OAB").grid(row=3, column=2, columnspan=self.lacols-4, sticky='nsew')

        for child in bottom.winfo_children(): child.grid_configure(padx=5, pady=5)    
        
        for child in laframe.winfo_children(): child.grid_configure(padx=5, pady=5)

    def update_text(self,message, timestamp=True):
        # Update the text in the log window
        if timestamp:
            now = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
            self.messages_entry.insert(tk.END, ': '.join((now,message)))
        else:
            self.messages_entry.insert(tk.END, message)
        self.messages_entry.yview(tk.END)  # Auto-scroll
        self.messages_entry.update_idletasks()
        return
        
    def reset(self):
        # reset all fields, clear log and plot
        self.set_entries()
        self.messages_entry.delete("1.0",tk.END)
        self.messages_entry.insert(tk.END, self.messages.get())
        self.plot_reset = True
        return


    def abort(self):
        # end all processes
        if self.thread_running.get():
            self.abort_value = True
            self.update_text("\n\n########################\nAborting the processes\nWaiting for the right moment...\n########################\n\n")
        else:
            self.update_text("\n\n########################\nNo process is running.\n########################\n\n")        
        return

        
    def quit_shiva(self):
        # quit the GUI
        self.abort()
        self.close_plot = True
        self.master.quit()
        return

    def nor_load_indir(self):
        # select normalisation input directory        
        loadbasedir = self.nor_indir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askdirectory(initialdir = loadbasedir,title = "Select input directory",mustexist=True)
        if fname:
            self.nor_indir.set(fname)
        else:
            fname = loadbasedir
            self.nor_indir.set(fname)
        self.nor_outdir.set(os.path.join(fname,self.outdir)) 
        self.lsd_indir.set(self.nor_outdir.get())
        self.lsd_outdir.set(self.nor_outdir.get())
        self.la_indir.set(self.nor_outdir.get())
        self.la_outdir.set(self.nor_outdir.get())
        return
        
    def nor_load_outdir(self):
        # select normalisation output directory
        loadbasedir = self.nor_outdir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askdirectory(initialdir = loadbasedir,title = "Select input directory")
        if fname:
            self.nor_outdir.set(fname)
        else:
            fname = loadbasedir
            self.nor_outdir.set(fname)
        self.lsd_indir.set(fname)
        self.la_indir.set(fname)
        self.lsd_outdir.set(fname)
        self.la_outdir.set(fname)
        return
    
    def change_cols(self, select_combobox):
        selection = self.option_instr.get()
        # self.option_values = ['UNDEF', 'ESPRESSO S1D', 'ESPRESSO S2D', 'GIANO-B MS1D', 'CARMENES', 'HARPS(N) S1D NEW DRS', 'HARPS(N) S2D NEW DRS'] 
        # self.units_values = ['angstroms','nm','micron']
        if selection == 'ESPRESSO S1D':
            self.wavecol.set(1)
            self.fluxcol.set(3)
            self.snrcol.set(0)
            self.errcol.set(4)
            self.echcol.set(0)
            self.spec_unit.set('a')
            self.units_default.set(self.units_values[0])
        elif selection == 'ESPRESSO S2D':
            self.wavecol.set(4)
            self.fluxcol.set(1)
            self.snrcol.set(0)
            self.errcol.set(2)
            self.echcol.set(0)
            self.spec_unit.set('a')
            self.units_default.set(self.units_values[0])
        elif selection == 'GIANO-B MS1D':
            self.wavecol.set(2)
            self.fluxcol.set(3)
            self.snrcol.set(4)
            self.errcol.set(0)
            self.echcol.set(1)
            self.spec_unit.set('n')
            self.units_default.set(self.units_values[1])
        elif selection == 'CARMENES':
            self.wavecol.set(4)
            self.fluxcol.set(1)
            self.snrcol.set(0)
            self.errcol.set(3)
            self.echcol.set(0)
            self.spec_unit.set('a')
            self.units_default.set(self.units_values[0])
        elif selection == 'HARPS(N) S1D NEW DRS':
            self.wavecol.set(1)
            self.fluxcol.set(3)
            self.snrcol.set(0)
            self.errcol.set(4)
            self.echcol.set(0)
            self.spec_unit.set('a')
            self.units_default.set(self.units_values[0])
        elif selection == 'HARPS(N) S2D NEW DRS':
            self.wavecol.set(4)
            self.fluxcol.set(1)
            self.snrcol.set(0)
            self.errcol.set(2)
            self.echcol.set(0)
            self.spec_unit.set('a')
            self.units_default.set(self.units_values[0])
        else:
            self.wavecol.set(1)
            self.fluxcol.set(2)
            self.snrcol.set(0)
            self.errcol.set(0)
            self.echcol.set(0)
            self.spec_unit.set('a')
            self.units_default.set(self.units_values[0])
        return
 
    def change_units(self, select_combobox):
        selection = self.units_default.get()
        self.spec_unit.set(selection[0])
        return
        
    def mask_change_units(self, select_combobox):
        selection = self.mask_units_default.get()
        self.mask_unit.set(selection[0])
        return
        
    def change_frame(self, select_combobox):
        selection = self.wave_frame.get()
        if selection == self.wave_values[0]:
            self.spec_vacuum.set(1)
        else:
            self.spec_vacuum.set(0)
        return
        
    def mask_change_frame(self, select_combobox):
        selection = self.mask_wave_frame.get()
        if selection == self.mask_wave_values[0]:
            self.mask_vacuum.set(1)
        else:
            self.mask_vacuum.set(0)
        return

    def option_la_change(self, select_combobox):
        selection = self.option_limits.get()
                
        if selection == self.option_la_values[0]:
            self.show_limitup.set(u'3 \u03c3')
            self.show_limitlow.set(u'3 \u03c3')
            self.usegauss.set(1)
            self.userot.set(0)
        elif selection == self.option_la_values[1]:
            self.show_limitup.set('vsini')
            self.show_limitlow.set('vsini')
            self.usegauss.set(0)
            self.userot.set(1)
        else:
            self.show_limitup.set(str(self.limitup.get()))
            self.show_limitlow.set(str(self.limitlow.get()))
            self.usegauss.set(0)
            self.userot.set(0)
        
        return

        
    def lsd_load_indir(self):
        # select LSD input directory
        loadbasedir = self.lsd_indir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askdirectory(initialdir = loadbasedir,title = "Select input directory",mustexist=True)
        if fname:
            self.lsd_indir.set(fname)
        else:
            fname = loadbasedir
            self.lsd_indir.set(fname)
        self.la_indir.set(fname)
        return
        
    def lsd_load_outdir(self):
        # select LSD output directory
        loadbasedir = self.lsd_outdir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askdirectory(initialdir = loadbasedir,title = "Select input directory")
        if fname:
            self.lsd_outdir.set(fname)
        else:
            fname = loadbasedir
            self.lsd_outdir.set(fname)
        self.la_indir.set(fname)
        self.la_outdir.set(fname)
        return
        
    def la_load_indir(self):
        # select Line Analysis input directory
        loadbasedir = self.la_indir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askdirectory(initialdir = loadbasedir,title = "Select input directory",mustexist=True)
        if fname:
            self.la_indir.set(fname)
        else:
            fname = loadbasedir
            self.la_indir.set(fname)
        return
        
    def la_load_outdir(self):
        # select Line Analysis output directory
        loadbasedir = self.la_outdir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askdirectory(initialdir = loadbasedir,title = "Select input directory")
        if fname:
            self.la_outdir.set(fname)
        else:
            fname = loadbasedir
            self.la_outdir.set(fname)
        return

    def nor_load_file(self):
        # load spectrum file
        loadbasedir = self.nor_indir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askopenfilename(initialdir = loadbasedir, title = "Select input file")
        if fname:
            self.nor_spec.set(fname)
            loadbasedir = os.path.split(fname)[0]
        else:
            fname = self.nor_spec.get()
            self.nor_spec.set(fname)

        self.nor_indir.set(loadbasedir)
        self.nor_outdir.set(os.path.join(loadbasedir,self.outdir)) 
        self.lsd_indir.set(self.nor_outdir.get())
        self.lsd_outdir.set(self.nor_outdir.get())
        self.la_indir.set(self.nor_outdir.get())
        self.la_outdir.set(self.nor_outdir.get())
        return

    def lsd_load_file(self):
        # load normalised spectrum file
        loadbasedir = self.lsd_indir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askopenfilename(initialdir = loadbasedir, title = "Select normalised file")
        if fname:
            self.lsd_spec.set(fname)
            self.lsd_indir.set(os.path.split(fname)[0])
        else:
            fname = self.lsd_spec.get()
            self.lsd_spec.set(fname)
            self.lsd_indir.set(loadbasedir)
        return
        
    def mask_load_file(self):
        # load mask file
        loadbasedir = self.nor_indir.get()
        checkdir = os.path.isdir(loadbasedir)
        if not checkdir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askopenfilename(initialdir = loadbasedir, title = "Select mask file")
        if fname:
            self.mask.set(fname)
        return
 

    def la_load_file(self):
        # load file for line analysis
        loadbasedir = os.path.isdir(self.la_indir.get())
        if not loadbasedir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        fname = filedialog.askopenfilename(initialdir = loadbasedir, title = "Select profile file")
        if fname:
            self.la_spec.set(fname)
            self.la_indir.set(os.path.split(fname)[0])

        return
        
    def save_log(self):
        # Save all the text in the log window to an ASCII file
        log = self.messages_entry.get("1.0",tk.END)
        date = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
        logname = '_'.join((date,'log_shiva'))
        loadbasedir = os.path.isdir(self.nor_indir.get())
        if not loadbasedir:
            loadbasedir = self.basedir
            if not os.path.isdir(loadbasedir):
                loadbasedir = os.getcwd()
        lname = filedialog.asksaveasfilename(initialdir = loadbasedir, initialfile = logname , defaultextension='.txt', confirmoverwrite=True, title = "Save log file")
        with open(lname, 'w') as l:
            l.write(log)
        return

    def plot_data(self):
        self.fig, self.ax = plt.subplots(figsize=(4,3), dpi=100)
        self.plot_logo()
        while True:
            time.sleep(0.01)
            if self.do_plot:
                self.update_figure(self.xvalues, self.yvalues, self.y_add, self.y_res, self.ymin, self.ymax, self.limits, self.logscale)
                self.do_plot = False
            if self.close_plot:
                return
            if self.plot_reset:
                self.plot_logo()
                self.plot_reset = False
        return


    def thread_plot(self):
        # thread running matplotlib
        self.th_plot = Thread(target=self.plot_data, daemon=True)
        self.th_plot.start()
        return
        
    def plot_logo(self):
        self.ax.clear()
        try:
            self.ax1.clear()
        except AttributeError:
            pass
        try:
            img = plt.imread(self.sp_logo)
            self.ax.imshow(img)
        except FileNotFoundError:
            img = np.zeros((40,30))
            self.ax.imshow(img, cmap='Greys')
        
        self.canvas = FigureCanvasTkAgg(self.fig, self.pltframe)
        self.canvas.get_tk_widget().grid(row=0, rowspan=self.pltrows, column=0, columnspan=self.pltcols, sticky='nsew') 
        self.fig.canvas.draw()
        return
            

    def update_figure(self, xvalues, yvalues, y_adds=None, y_res=None, ymin=None, ymax=None, limits=(None,None), logscale=False):
        # Update the figure in the plot window
        plt.clf()
        plt.close(self.fig)
       

        if y_res is None:
            self.fig, self.ax = plt.subplots(figsize=(4,3), dpi=100)
        else:
            self.fig, [self.ax, self.ax1] = plt.subplots(nrows=2, ncols=1, figsize=(4,3), dpi=100, gridspec_kw={'height_ratios': [3, 1]})

        if ymin is not None:
            self.ax.set_ylim(bottom=ymin)
        if ymax is not None:
            self.ax.set_ylim(top=ymax)
                        
        self.ax.plot(xvalues,yvalues)
        if y_adds is not None:
            for y_add in y_adds:
                self.ax1.plot(xvalues,y_add)        
        for limit in limits:
            if limit is not None:
                self.ax.axvline(limit)

        if logscale:
            plt.gca().set_yscale('log')
        else:
            plt.gca().set_yscale('linear')
        self.canvas = FigureCanvasTkAgg(self.fig, self.pltframe)
        self.canvas.get_tk_widget().grid(row=0, rowspan=self.pltrows, column=0, columnspan=self.pltcols, sticky='nsew')
        self.fig.canvas.draw()
        #print(2)
        return

    def find_file(self,filepattern,fileplace=None):
        # Find the input file
        if os.path.isfile(filepattern):
            return [filepattern]
        else:
            wildcards = ''.join(('*',filepattern))
            if fileplace is None:
                fileplace = self.basedir
            specs = glob.glob(os.path.join(fileplace,filepattern))
            specs = sorted(specs)
            return specs
        
    def save_fits(self,data,header,oldname,outdir,suffix):
        # Save the output as a FITS BinTable
        columns = []
        for key in data:
            if isinstance(data[key], np.ndarray):
                columns.append(fits.Column(name=key.upper(), format='D', array=data[key]))
            else:
                columns.append(fits.Column(name=key.upper(), format='D', array=np.ones(2)*data[key]))
        tbhdu = fits.BinTableHDU.from_columns(columns)
        prihdu = fits.PrimaryHDU(data=None, header=header)
        hdulist = fits.HDUList([prihdu, tbhdu])

        if not os.path.isdir(outdir):
            os.mkdir(outdir)
            
        basename = os.path.basename(oldname)
        basename = os.path.splitext(basename)[0]
        newname = ''.join((basename,suffix))
        newname = os.path.join(outdir,newname)
        hdulist.writeto(newname,overwrite=True,output_verify='ignore')
        
        return
        
    def insert_key(self, hea, keyword, value, comment):
        # check if keyword value is allowed (e.g., not NaN)
        value = np.nan_to_num(value, nan=0.0, posinf=0.0, neginf=0.0)
        hea[keyword] = (value, comment)
        return hea

    def update_fits(self,fitsname,hea):
        # update header FITS
        with fits.open(fitsname, mode='update') as hdu:
            for entry in hea:
                try:
                    hdu[0].header[entry] = hea[entry]
                except ValueError:
                    pass
            hdu.flush()
        return

        
    def sfx_lsd(self):
        self.lsd_output.set('_lsd.fits')    # suffix of the output LSD profiles
        self.do_ccf.set(0)
        self.la_spec.set('*_lsd.fits')      # input (normalised) LSD profile OR pattern
        return

    def sfx_ccf(self):
        self.lsd_output.set('_ccf.fits')    # suffix of the output CCF profiles
        self.do_lsd.set(0)
        self.la_spec.set('*_ccf.fits')      # input (normalised) CCF profile OR pattern
        return
    
    def thread_finished(self):
        self.thread_running.set(0)
        return
        
    def normalise(self):
        # Normalise spectra
        
        self.update_text("\n#######################\n", timestamp=False)
        self.update_text("   START normalising the spectra\n", timestamp=False)
        self.update_text("#######################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Units: {self.spec_unit.get()}\n", timestamp=False)
        self.update_text(f"Polynomial degree: {self.degree.get()}\n", timestamp=False)
        self.update_text(f"Echelle: {self.echcol.get()} (if 0, subsets are used)\n", timestamp=False)
        self.update_text(f"Subsets: {self.n_ord.get()} (if 0, whole spectrum is used)\n", timestamp=False)
        self.update_text(f"Refine: {self.refine.get()}\n\n", timestamp=False)
        spectra = self.find_file(self.nor_spec.get(), self.nor_indir.get())
        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n#######################\n", timestamp=False)
            self.update_text("   END normalising the spectra\n", timestamp=False)
            self.update_text("#######################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} spectra. Start normalising.\n")
        # iter on spectra, normalise, save output, plot spectra, update variables
        for n, spec in enumerate(spectra):
            if self.abort_value:
                self.update_text("Normalisation aborted.\n")
                self.abort_value = False
                break
            self.update_text(f"{n+1} of {n_spectra} spectra: working on spectrum {os.path.basename(spec)}.\n")
            spectrum = pa.read_spectrum(spec, unit=self.spec_unit.get(), wavecol=self.wavecol.get(), \
                 fluxcol=self.fluxcol.get(), snrcol=self.snrcol.get(), echcol=self.echcol.get(), errcol=self.errcol.get(), vacuum=bool(self.spec_vacuum.get()))
            nor = pa.norm_spectrum(spectrum['wave'], spectrum['flux'], spectrum['snr'], spectrum['echelle'],\
                deg=self.degree.get(), n_ord=self.n_ord.get(), refine=self.refine.get(), output=False)
            hea = spectrum['header']
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_nor[0], self.degree.get(),self.key_nor[1])
            if self.n_ord.get():
                hea = self.insert_key(hea , self.key_sub[0] , self.n_ord.get(),self.key_sub[1])
            else:
                try:
                    subsets = max(spectrum['echelle']) - min(spectrum['echelle']) +1
                except TypeError:
                    subsets = self.n_ord.get()
                hea = self.insert_key(hea, self.key_sub[0], subsets,self.key_sub[1])
            hea = self.insert_key(hea, self.key_refine[0], bool(self.refine.get()), self.key_refine[1])
            hea = self.insert_key(hea, self.key_swave[0], self.spec_unit.get(), self.key_swave[1])
            self.save_fits(nor,hea,spec,self.nor_outdir.get(),self.nor_output.get())
            
            if self.norplot.get():
                self.xvalues = nor['wave']
                self.yvalues = nor['nflux']
                self.y_add = None
                self.y_res = None
                self.ymin = 0
                self.ymax = 1.5 
                self.limits = (None,None)
                self.logscale = False
                self.do_plot = True
            self.update_text(f"Spectrum {n+1} normalised.\n")

 
        self.update_text("\n#######################\n", timestamp=False)
        self.update_text("   END normalising the spectra\n", timestamp=False)
        self.update_text("#######################\n\n", timestamp=False)
        #self.thread_finished()
        return

    def call_normalise(self):
        try:
            self.normalise()
        except:
            self.update_text("\n\n########################\nAborting the process.\nAre there any spectra with the same extension but different format?\n########################\n\n")
        self.thread_finished()        
        return

    def thread_normalise(self):
        # this thread call the "normalise" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_normalise = Thread(target=self.call_normalise, daemon=True)
            self.thread_running.set(1)
            self.th_normalise.start()
        return

    def do_profile(self):
        # Create profile (LSD or CCF)
        self.update_text("\n###########################\n", timestamp=False)
        if self.do_ccf.get():
            self.update_text("  START computing CCF profiles \n", timestamp=False)
            if self.ccfweight.get():
                self.update_text("  Using the S/N as weigths \n", timestamp=False)
            root_key = self.base_key['ccf']
        else:
            self.update_text("  START computing LSD profiles \n", timestamp=False)
            root_key = self.base_key['lsd']
        self.update_text("###########################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Spectral units: {self.lsd_spec_unit.get()}\n", timestamp=False)
        self.update_text(f"Absorption mask: {self.mask_invert.get()}\n", timestamp=False)
        self.update_text(f"Mask units: {self.mask_unit.get()}\n", timestamp=False)
        self.update_text(f"Mask select line from depths: ({self.mask_dlow.get()},{self.mask_dup.get()})\n", timestamp=False)
        self.update_text(f"Mask select line from wavelength: ({self.mask_wmin.get()},{self.mask_wmax.get()})\n", timestamp=False)
        self.update_text(f"Exclude Balmer regions: {bool(self.mask_balmer.get())}\n", timestamp=False)
        self.update_text(f"Exclude telluric regions: {bool(self.mask_tell.get())}\n", timestamp=False)
        self.update_text(f"Remove cosmics: {bool(self.cosmic.get())}\n", timestamp=False)
        self.update_text(f"Clean profile: {bool(self.clean.get())}\n", timestamp=False)
        self.update_text(f"RV range and step: ({self.rvmin.get()},{self.rvmax.get()}) {self.rvstep.get()}\n\n", timestamp=False)
        spectra = self.find_file(self.lsd_spec.get(), self.lsd_indir.get())

        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n###########################\n", timestamp=False)
            if self.do_ccf.get():
               self.update_text("  END computing CCF profiles \n", timestamp=False)
            else:
               self.update_text("  END computing LSD profiles \n", timestamp=False)
            self.update_text("###########################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} spectra.\n")
        
        # Find the mask file and read it
        maskfile = self.find_file(self.mask.get())[0]
        try:
            mask = pa.read_mask(maskfile, unit = self.mask_unit.get(), ele = self.mask_els.get(), \
               no_ele = self.mask_noels.get(), depths=(self.mask_dlow.get(),self.mask_dup.get()), \
               balmer = bool(self.mask_balmer.get()), tellurics = bool(self.mask_tell.get()), \
               wmin = self.mask_wmin.get(), wmax = self.mask_wmax.get(), \
               invert = bool(self.mask_invert.get()), vacuum=bool(self.mask_vacuum.get()))
        except IsADirectoryError:
               self.update_text(f"No mask was defined. Aborted.\n")
               self.update_text("\n###########################\n", timestamp=False)
               if self.do_ccf.get():
                   self.update_text("  END computing CCF profiles \n", timestamp=False)
               else:
                   self.update_text("  END computing LSD profiles \n", timestamp=False)
               self.update_text("###########################\n\n", timestamp=False)
               return
        self.update_text(f"Read mask {os.path.basename(maskfile)}.\n") 
        
        # Iterate on normalised spectra: read the data, compute LSD, plot and save
        for n, spec in enumerate(spectra):
            if self.abort_value:
                self.update_text("Profile computation aborted.\n")
                self.abort_value = False
                break
            self.update_text(f"Working on spectrum {n+1} of {n_spectra}.\n")
            spectrum = pa.read_spectrum(spec, unit=self.lsd_spec_unit.get(), wavecol=self.lsd_wavecol.get(), \
                 fluxcol=self.lsd_fluxcol.get(), snrcol=self.lsd_snrcol.get(), nfluxcol=self.lsd_nfluxcol.get())
            try:
                with fits.open(spec) as hdu:
                    hea = hdu[0].header
                for jd_key in self.input_jd:
                    try:
                        jd = hea[jd_key]
                        ori_jd_key = jd_key
                        break
                    except KeyError:
                        jd = False
                        ori_jd_key = 'NONE'
            except OSError:
                jd = False
                ori_jd_key = 'NONE'
                
            
            ### DEFINE LSD OR CCF
            if self.do_ccf.get():
                profile = pa.compute_ccf(spectrum, mask, vrange=(self.rvmin.get(),self.rvmax.get()), step=self.rvstep.get(), mask_spectrum = bool(self.mask_spectrum.get()), cosmic=bool(self.cosmic.get()), clean=bool(self.clean.get()), weights=bool(self.ccfweight.get()), verbose=False, output=False)
            else:
                profile = pa.compute_lsd(spectrum, mask, vrange=(self.rvmin.get(),self.rvmax.get()), step=self.rvstep.get(), cosmic=bool(self.cosmic.get()), clean=bool(self.clean.get()), verbose=False, output=False)

            hea = fits.PrimaryHDU().header
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_lsdin[0])), os.path.basename(spec), self.key_lsdin[1])
            hea = self.insert_key(hea, self.key_jd[0], jd, ''.join((self.key_jd[1],ori_jd_key)))
            hea = self.insert_key(hea, ' '.join((root_key,self.key_mask[0])), os.path.basename(maskfile), self.key_mask[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_mask_invert[0])), bool(self.mask_invert.get()), self.key_mask_invert[1])            
            hea = self.insert_key(hea, ' '.join((root_key,self.key_mask_spectrum[0])), bool(self.mask_spectrum.get()), self.key_mask_spectrum[1])            
            
            hea = self.insert_key(hea, ' '.join((root_key,self.key_els[0])), self.mask_els.get(), self.key_els[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_noels[0])), self.mask_noels.get(), self.key_noels[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_dlow[0])), self.mask_dlow.get(), self.key_dlow[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_dup[0])), self.mask_dup.get(), self.key_dup[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_wmin[0])), self.mask_wmin.get(), self.key_wmin[1])   
            hea = self.insert_key(hea, ' '.join((root_key,self.key_wmax[0])), self.mask_wmax.get(), self.key_wmax[1])  
            hea = self.insert_key(hea, ' '.join((root_key,self.key_balmer[0])), bool(self.mask_balmer.get()), self.key_balmer[1])   
            hea = self.insert_key(hea, ' '.join((root_key,self.key_tell[0])), bool(self.mask_tell.get()), self.key_tell[1])         
            hea = self.insert_key(hea, ' '.join((root_key,self.key_cosmic[0])), bool(self.cosmic.get()), self.key_cosmic[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_clean[0])), bool(self.clean.get()), self.key_clean[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_rvmin[0])), self.rvmin.get(), self.key_rvmin[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_rvmax[0])), self.rvmax.get(), self.key_rvmax[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_rvstep[0])), self.rvstep.get(), self.key_rvstep[1])
            hea = self.insert_key(hea, self.key_ccfweight[0], bool(self.ccfweight.get()), self.key_ccfweight[1])
            self.save_fits(profile,hea,spec,self.lsd_outdir.get(),self.lsd_output.get())
            
            if self.lsdplot.get():
                self.xvalues = profile['rv_range']
                self.yvalues = profile['profile']
                self.y_add = None
                self.y_red = None
                self.ymin = None
                self.ymax = None
                self.limits = (None,None)
                self.logscale = False
                self.do_plot = True
            
                #plt.close(self.fig)
                #self.update_figure(profile['rv_range'], profile['profile'])

            self.update_text(f"Profile {n+1} computed.\n")
        self.limitlow.set(self.rvmin.get()+20)
        self.limitup.set(self.rvmax.get()-20)
        self.update_text("\n###########################\n", timestamp=False)
        if self.do_ccf.get():
            self.update_text("  END computing CCF profiles \n", timestamp=False)
        else:
            self.update_text("  END computing LSD profiles \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        #self.thread_finished()

        return

    def call_profile(self):
        try:
            self.do_profile()
        except:
            self.update_text("\n\n########################\nAborting the process.\nAre there any spectra with the same extension but different format?\n########################\nIs the mask/model covering at least part of the spectral range of the spectra?\n\n")
        self.thread_finished()        
        return

    def thread_do_profile(self):
        # this thread call the "do_profile" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_do_profile = Thread(target=self.call_profile, daemon=True)
            self.thread_running.set(1)
            self.th_do_profile.start()
        return

    def extract_line(self):
        # Extract a single line
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  START extracting single line \n", timestamp=False)
        root_key = self.base_key['line']

        self.update_text("###########################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Spectral units: {self.lsd_spec_unit.get()}\n", timestamp=False)
        self.update_text(f"Central wavelength of line: {self.ext_wave.get()}\n", timestamp=False)
        self.update_text(f"RV range and step: ({self.rvmin.get()},{self.rvmax.get()}) {self.rvstep.get()}\n\n", timestamp=False)

        spectra = self.find_file(self.lsd_spec.get(), self.lsd_indir.get())

        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n###########################\n", timestamp=False)
            self.update_text("  END extracting single line \n", timestamp=False)
            self.update_text("###########################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} spectra.\n")

        
        # Iterate on normalised spectra: read the data, extract the line, plot and save
        for n, spec in enumerate(spectra):
            if self.abort_value:
                self.update_text("Line extraction aborted.\n")
                self.abort_value = False
                break
            self.update_text(f"\nWorking on spectrum {n+1} of {n_spectra}.\n")
            spectrum = pa.read_spectrum(spec, unit=self.lsd_spec_unit.get(), wavecol=self.lsd_wavecol.get(), \
                 fluxcol=self.lsd_fluxcol.get(), snrcol=self.lsd_snrcol.get(), nfluxcol=self.lsd_nfluxcol.get())
            try:
                with fits.open(spec) as hdu:
                    hea = hdu[0].header
                for jd_key in self.input_jd:
                    try:
                        jd = hea[jd_key]
                        ori_jd_key = jd_key
                        break
                    except KeyError:
                        jd = False
                        ori_jd_key = 'NONE'
            except OSError:
                jd = False
                ori_jd_key = 'NONE'
                
            
            ### Extract the line
            profile = pa.extract_line(spectrum, unit=self.lsd_spec_unit.get(), w0=self.ext_wave.get(), vrange=(self.rvmin.get(),self.rvmax.get()), step=self.rvstep.get(), verbose=False, output=False)
            
            hea = fits.PrimaryHDU().header
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_lsdin[0])), os.path.basename(spec), self.key_lsdin[1])
            hea = self.insert_key(hea, self.key_jd[0], jd, ''.join((self.key_jd[1],ori_jd_key)))
            hea = self.insert_key(hea, ' '.join((root_key,self.key_rvmin[0])), self.rvmin.get(), self.key_rvmin[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_rvmax[0])), self.rvmax.get(), self.key_rvmax[1])
            hea = self.insert_key(hea, ' '.join((root_key,self.key_rvstep[0])), self.rvstep.get(), self.key_rvstep[1])
            self.save_fits(profile,hea,spec,self.lsd_outdir.get(),self.ext_output.get())
            
            if self.extplot.get():
                self.xvalues = profile['rv_range']
                self.yvalues = profile['profile']
                self.y_add = None
                self.y_red = None
                self.ymin = None
                self.ymax = None
                self.limits = (None,None)
                self.logscale = False
                self.do_plot = True

            self.update_text(f"Line {n+1} extracted.\n")
        self.limitlow.set(self.rvmin.get()+20)
        self.limitup.set(self.rvmax.get()-20)
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  END extracting single line \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        #self.thread_finished()
        
        return    

    def call_extract(self):
        try:
            self.extract_line()
        except:
            self.update_text("\n\n########################\nAborting the process.\n- Are there any spectra with the same extension but different format?\n- Maybe the extraction window falls on overlapping orders.\n  If the echelle orders are defined in the data, try to restrict the window size.\n  Otherwise, re-normalise the spectra with subsets=0 and degree=0 (tab 1)\n  You can then normalise the profiles in tab 3a\n########################\n\n")
        self.thread_finished()        
        return
    
    def thread_extract_line(self):
        # this thread call the "extract_profile" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_extract_line = Thread(target=self.call_extract, daemon=True)
            self.thread_running.set(1)
            self.th_extract_line.start()
        return


    def norm_profile(self):
        # Normalise profiles
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  START normalising profiles \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Define continuum: ({self.limitlow.get()},{self.limitup.get()})\n", timestamp=False)

        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n###########################\n", timestamp=False)
            self.update_text("  END normalising profiles \n", timestamp=False)
            self.update_text("###########################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} profiles.\n")
        
        if self.std.get():
            stdbasedir = self.la_outdir.get()
            stdbasename = os.path.join(stdbasedir,'line_mean_std')
            pfns, stds = pa.norm_profile(spectra, rvcol=1, prfcol=2, errcol=3, sfx=False, std=stdbasename, limits=(self.limitlow.get(),self.limitup.get()))
        else:
            pfns, stds = pa.norm_profile(spectra, rvcol=1, prfcol=2, errcol=3, sfx=False, std=False, limits=(self.limitlow.get(),self.limitup.get()))
        
        for n,pfn in enumerate(pfns):
            if self.abort_value:
                self.update_text("Profile normalisation aborted.\n")
                self.abort_value = False
                break
            hea = pfn.pop('header')
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spectra[n]), self.key_lainp[1])
            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lnor[0], True, self.key_lnor[1])
            hea = self.insert_key(hea, self.key_lnrvmin[0], self.limitlow.get(), self.key_lnrvmin[1])
            hea = self.insert_key(hea, self.key_lnrvmax[0], self.limitup.get(), self.key_lnrvmax[1])
            self.save_fits(pfn,hea,spectra[n],self.la_outdir.get(),self.norprf_output.get())
            
            if self.norprfplot.get():
                self.xvalues = pfn['rv_range']
                self.yvalues = pfn['profile']
                self.y_add = [pfn['nprofile']]
                self.y_res = None
                self.ymin = None
                self.ymax = None
                self.limits = (None,None)
                self.logscale = False
                self.do_plot = True

            self.update_text(f"Profile {n+1} computed.\n")            
        if self.std.get():
            self.xvalues = stds['rv_mean']
            self.yvalues = stds['ccf_mean']
            self.y_add = None
            self.y_res = stds['std_dev']
            self.ymin = None
            self.ymax = None
            self.limits = (None,None)
            self.logscale = False
            self.do_plot = True

        self.la_spec.set(''.join(('*',self.norprf_output.get())))

        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  END normalising profiles \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        #self.thread_finished()
        return

    def call_norprofile(self):
        try:
            self.norm_profile()
        except:
            self.update_text("\n\n########################\nAborting the process.\nDo the profiles have the same RV range?\n########################\n\n")
        self.thread_finished()        
        return

    def thread_norm_profile(self):
        # this thread call the "norm_profile" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_norm_profile = Thread(target=self.call_norprofile, daemon=True)
            self.thread_running.set(1)
            try:
                self.th_norm_profile.start()
            except:
                print('Exception')
                self.thread_running.set(0)
        return

         
    def fit_profile(self):
        # Fit profiles
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  START fitting profiles \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Use Gaussian fitting: {self.fitgauss.get()},\n", timestamp=False)
        self.update_text(f"Use Lorentzian fitting: {self.fitlorentz.get()},\n", timestamp=False)
        self.update_text(f"Use Voigt fitting: {self.fitvoigt.get()},\n", timestamp=False)
        self.update_text(f"Use rotational fitting: {self.fitrot.get()},\n", timestamp=False)
        self.update_text(f"Linear Limb Darkening: {self.ld.get()},\n", timestamp=False)
        self.update_text(f"RV guess value: {self.rv0.get()},\n", timestamp=False)
        self.update_text(f"Width guess value: {self.width.get()},\n\n", timestamp=False)

        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n###########################\n", timestamp=False)
            self.update_text("  END fitting profiles \n", timestamp=False)
            self.update_text("###########################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} profiles.\n")
        
        for n, spec in enumerate(spectra):
            if self.abort_value:
                self.update_text("Profile fitting aborted.\n")
                self.abort_value = False
                break
            self.update_text(f"{n+1} of {n_spectra} profiles: working on profile {os.path.basename(spec)}.\n")
            with fits.open(spec) as hdu:
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            vmax = np.amax(vrad)
            vmin = np.amin(vrad)
            if np.logical_or(self.rv0.get()<= vmin, self.rv0.get()>= vmax):
                self.rv0.set(int((vmax+vmin)/2))
                self.update_text("********* WARNING *********\n",timestamp=False)
                self.update_text(f"RV guess value outside the RV range. It was changed to {self.rv0.get()}.\n")
                self.update_text("********* END WARNING *********\n\n",timestamp=False)
            flux = data.field(1)
            if self.fit_errs.get():
                errs = data.field(2)
            else:
                errs = 0
            fit_result = pa.fit_profile(vrad, flux, errs=errs, gauss=bool(self.fitgauss.get()), rot=bool(self.fitrot.get()), lorentz=bool(self.fitlorentz.get()), voigt=bool(self.fitvoigt.get()), rv0=self.rv0.get(), width=self.width.get(), ld=self.ld.get())
            if self.fitgauss.get():
                self.update_text(f"Gaussian fit:\n")
                self.update_text(f"RV = {np.round(fit_result['gaussian']['rv'],4)} +/- {np.round(fit_result['gaussian']['e_rv'],4)}  km/s\n", timestamp=False)
                self.update_text(f"FWHM = {np.round(fit_result['gaussian']['fwhm'],4)} +/- {np.round(fit_result['gaussian']['e_fwhm'],4)}  km/s\n", timestamp=False)
                self.update_text(f"EW = {np.round(fit_result['gaussian']['EW'],4)} +/- {np.round(fit_result['gaussian']['e_EW'],4)}  km/s\n\n", timestamp=False)
            if self.fitlorentz.get():
                self.update_text(f"Lorenztian fit:\n")
                self.update_text(f"RV = {np.round(fit_result['lorentzian']['rv'],4)} +/- {np.round(fit_result['lorentzian']['e_rv'],4)}  km/s\n", timestamp=False)
                self.update_text(f"FWHM = {np.round(fit_result['lorentzian']['fwhm'],4)} +/- {np.round(fit_result['lorentzian']['e_fwhm'],4)}  km/s\n", timestamp=False)
                self.update_text(f"EW = {np.round(fit_result['lorentzian']['EW'],4)} +/- {np.round(fit_result['lorentzian']['e_EW'],4)}  km/s\n\n", timestamp=False)
            if self.fitvoigt.get():
                self.update_text(f"Voigt fit:\n")
                self.update_text(f"RV = {np.round(fit_result['voigt']['rv'],4)} +/- {np.round(fit_result['voigt']['e_rv'],4)}  km/s\n", timestamp=False)
                self.update_text(f"FWHM = {np.round(fit_result['voigt']['fwhm'],4)} +/- {np.round(fit_result['voigt']['e_fwhm'],4)}  km/s\n", timestamp=False)
                self.update_text(f"EW = {np.round(fit_result['voigt']['EW'],4)} +/- {np.round(fit_result['voigt']['e_EW'],4)}  km/s\n\n", timestamp=False)
            if self.fitrot.get():
                self.update_text(f"Rotational fit:\n")
                self.update_text(f"RV = {np.round(fit_result['rotational']['rv'],4)} +/- {np.round(fit_result['rotational']['e_rv'],4)}  km/s\n", timestamp=False)
                self.update_text(f"Vsini = {np.round(fit_result['rotational']['vsini'],4)} +/- {np.round(fit_result['rotational']['e_vsini'],4)}  km/s\n", timestamp=False)
                self.update_text(f"EW = {np.round(fit_result['rotational']['EW'],4)} +/- {np.round(fit_result['rotational']['e_EW'],4)}  km/s\n", timestamp=False)

            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])      
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spec), self.key_lainp[1])
            hea = self.insert_key(hea, self.key_fitgauss[0], bool(self.fitgauss.get()), self.key_fitgauss[1])
            hea = self.insert_key(hea, self.key_fitrot[0], bool(self.fitrot.get()), self.key_fitrot[1])
            hea = self.insert_key(hea, self.key_fitld[0], self.ld.get(), self.key_fitld[1])
            hea = self.insert_key(hea, self.key_rvguess[0], self.rv0.get(), self.key_rvguess[1])
            hea = self.insert_key(hea, self.key_wguess[0], self.width.get(), self.key_wguess[1])
            hea = self.insert_key(hea, self.key_gaussrv[0], fit_result['gaussian']['rv'], self.key_gaussrv[1])
            hea = self.insert_key(hea, self.key_gaussrverr[0], fit_result['gaussian']['e_rv'], self.key_gaussrverr[1])
            hea = self.insert_key(hea, self.key_gaussfwhm[0], fit_result['gaussian']['fwhm'], self.key_gaussfwhm[0])
            hea = self.insert_key(hea, self.key_gaussfwhmerr[0], fit_result['gaussian']['e_fwhm'], self.key_gaussfwhmerr[1])
            hea = self.insert_key(hea, self.key_gaussew[0], fit_result['gaussian']['EW'], self.key_gaussew[1])
            hea = self.insert_key(hea, self.key_gaussewerr[0], fit_result['gaussian']['e_EW'], self.key_gaussewerr[1])
            hea = self.insert_key(hea, self.key_lorentzrv[0], fit_result['lorentzian']['rv'], self.key_lorentzrv[1])
            hea = self.insert_key(hea, self.key_lorentzrverr[0], fit_result['lorentzian']['e_rv'], self.key_lorentzrverr[1])
            hea = self.insert_key(hea, self.key_lorentzfwhm[0], fit_result['lorentzian']['fwhm'], self.key_lorentzfwhm[0])
            hea = self.insert_key(hea, self.key_lorentzfwhmerr[0], fit_result['lorentzian']['e_fwhm'], self.key_lorentzfwhmerr[1])
            hea = self.insert_key(hea, self.key_lorentzew[0], fit_result['lorentzian']['EW'], self.key_lorentzew[1])
            hea = self.insert_key(hea, self.key_lorentzewerr[0], fit_result['lorentzian']['e_EW'], self.key_lorentzewerr[1])
            hea = self.insert_key(hea, self.key_voigtrv[0], fit_result['voigt']['rv'], self.key_voigtrv[1])
            hea = self.insert_key(hea, self.key_voigtrverr[0], fit_result['voigt']['e_rv'], self.key_voigtrverr[1])
            hea = self.insert_key(hea, self.key_voigtfwhm[0], fit_result['voigt']['fwhm'], self.key_voigtfwhm[0])
            hea = self.insert_key(hea, self.key_voigtfwhmerr[0], fit_result['voigt']['e_fwhm'], self.key_voigtfwhmerr[1])
            hea = self.insert_key(hea, self.key_voigtew[0], fit_result['voigt']['EW'], self.key_voigtew[1])
            hea = self.insert_key(hea, self.key_voigtewerr[0], fit_result['voigt']['e_EW'], self.key_voigtewerr[1])
            hea = self.insert_key(hea, self.key_rotrv[0], fit_result['rotational']['rv'], self.key_rotrv[1])
            hea = self.insert_key(hea, self.key_rotrverr[0], fit_result['rotational']['e_rv'], self.key_rotrverr[1])
            hea = self.insert_key(hea, self.key_rotvsini[0], fit_result['rotational']['vsini'], self.key_rotvsini[1])
            hea = self.insert_key(hea, self.key_rotvsinierr[0], fit_result['rotational']['e_vsini'], self.key_rotvsinierr[1])
            hea = self.insert_key(hea, self.key_rotew[0], fit_result['rotational']['EW'], self.key_rotew[1])
            hea = self.insert_key(hea, self.key_rotewerr[0], fit_result['rotational']['e_EW'], self.key_rotewerr[1])
                        
            res_dict= {'rv_range' : vrad, 'profile' : flux, 'fitgauss' : fit_result['gaussian']['profile'], 'fitlorentz' : fit_result['lorentzian']['profile'], 'fitvoigt' : fit_result['voigt']['profile'], 'fitrot' : fit_result['rotational']['profile']}
            self.save_fits(res_dict,hea,spec,self.la_outdir.get(),self.fit_output.get())
            self.update_fits(spec,hea)
            if self.fitplot.get():
                self.xvalues = vrad
                self.yvalues = flux
                self.y_add = [fit_result['gaussian']['profile'], fit_result['lorentzian']['profile'], fit_result['voigt']['profile'], fit_result['rotational']['profile']]
                self.y_res = None
                self.ymin = None
                self.ymax = None
                self.limits = (None,None)
                self.logscale = False
                self.do_plot = True

        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  END fitting profiles \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        #self.thread_finished()
        return

    def call_fit(self):
        try:
            self.fit_profile()
        except:
            self.update_text("\n\n########################\nAborting the process.\nAre there any profiles with the same extension but different format?\n########################\n\n")
        self.thread_finished()        
        return

    def thread_fit_profile(self):
        # this thread call the "fit_profile" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_fit_profile = Thread(target=self.call_fit, daemon=True)
            self.thread_running.set(1)
            self.th_fit_profile.start()
        return
         
    def moments(self):
        # Compute moments
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  START computing the moments \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Use Gaussian fitting: {self.usegauss.get()},\n", timestamp=False)
        self.update_text(f"Use rotational fitting: {self.userot.get()},\n", timestamp=False)
        self.update_text(f"Manual limits: ({self.limitlow.get()},{self.limitup.get()})  km/s\n", timestamp=False)

        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n###########################\n", timestamp=False)
            self.update_text("  END computing the moments \n", timestamp=False)
            self.update_text("###########################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} profiles.\n")
        
        for n, spec in enumerate(spectra):
            if self.abort_value:
                self.update_text("Moments computation aborted.\n")
                self.abort_value = False
                break
            self.update_text(f"{n+1} of {n_spectra} profiles: working on profile {os.path.basename(spec)}.\n")
            with fits.open(spec) as hdu:
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            flux = data.field(1)
            errs = data.field(2)
            if self.usegauss.get():
                x0 = hea[self.key_gaussrv[0]]
                fwhm = hea[self.key_gaussfwhm[0]]
                sigma = fwhm/np.sqrt(8*np.log(2))
                limits = (x0-3*sigma, x0+3*sigma)
            elif self.userot.get():
                x0 = hea[self.key_rotrv[0]]
                vsini = hea[self.key_rotvsini[0]]
                limits = (x0-vsini, x0+vsini)
            else:
                limits = (self.limitlow.get(), self.limitup.get())
            mom = pa.moments(vrad, flux, errs=errs, limits=limits, normalise=True)
            self.update_text(f"Moments:\n")
            self.update_text(f"m0 (EW) = {round(mom['m0'],4)} +/- {round(mom['e_m0'],4)}  km/s\n", timestamp=False)
            self.update_text(f"m1 (RV) = {round(mom['m1'],4)} +/- {round(mom['e_m1'],4)}  km/s\n", timestamp=False)
            self.update_text(f"m2 (Variance) = {round(mom['m2'],4)} +/- {round(mom['e_m2'],4)}\n", timestamp=False)
            self.update_text(f"m3 (Seed for Skewness) = {round(mom['m3'],4)} +/- {round(mom['e_m3'],4)}\n", timestamp=False)
            self.update_text(f"m4 (Seed for Kurtosis) = {round(mom['m4'],4)} +/- {round(mom['e_m4'],4)}\n", timestamp=False)
            self.update_text(f"FWHM (from m2) = {round(mom['fwhm'],4)} +/- {round(mom['e_fwhm'],4)}\n", timestamp=False)
            self.update_text(f"Skewness (from m3) = {round(mom['skewness'],4)} +/- {round(mom['e_skewness'],4)}\n", timestamp=False)
            self.update_text(f"Kurtosis (from m4) = {round(mom['kurtosis'],4)} +/- {round(mom['e_kurtosis'],4)}\n", timestamp=False)

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
            hea = self.insert_key(hea, self.key_momskew[0], mom['skewness'], self.key_momskew[1])
            hea = self.insert_key(hea, self.key_momskewerr[0], mom['e_skewness'], self.key_momskewerr[1])
            hea = self.insert_key(hea, self.key_momkurt[0], mom['kurtosis'], self.key_momkurt[1])
            hea = self.insert_key(hea, self.key_momkurterr[0], mom['e_kurtosis'], self.key_momkurterr[1])            
            
            self.save_fits(mom,hea,spec,self.la_outdir.get(),self.mom_output.get())
            self.update_fits(spec,hea)
                        
            if self.momplot.get():
                self.xvalues = vrad
                self.yvalues = flux
                self.y_add = None
                self.y_res = None
                self.ymin = None
                self.ymax = None
                self.limits = limits
                self.logscale = False
                self.do_plot = True

        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  END computing the moments \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        #self.thread_finished()
        return

    def call_moments(self):
        try:
            self.moments()
        except:
            self.update_text("\n\n########################\nAborting the process.\nAre there any profiles with the same extension but different format?\n########################\n\n")
        self.thread_finished()        
        return

    def thread_moments(self):
        # this thread call the "moments" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_moments = Thread(target=self.call_moments, daemon=True)
            self.thread_running.set(1)
            self.th_moments.start()
        return
         
    def bisector(self):
        # Compute bisector
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  START computing the bisector \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Use Gaussian fitting: {self.usegauss.get()},\n", timestamp=False)
        self.update_text(f"Use rotational fitting: {self.userot.get()},\n", timestamp=False)
        self.update_text(f"Manual limits: ({self.limitlow.get()},{self.limitup.get()})  km/s\n", timestamp=False)

        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n###########################\n", timestamp=False)
            self.update_text("  END computing the bisector \n", timestamp=False)
            self.update_text("###########################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} profiles.\n")
        
        for n, spec in enumerate(spectra):
            if self.abort_value:
                self.update_text("Bisector computation aborted.\n")
                self.abort_value = False
                break
            self.update_text(f"{n+1} of {n_spectra} profiles: working on profile {os.path.basename(spec)}.\n")
            with fits.open(spec) as hdu:
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            flux = data.field(1)
            if self.usegauss.get():
                x0 = hea[self.key_gaussrv[0]]
                fwhm = hea[self.key_gaussfwhm[0]]
                sigma = fwhm/np.sqrt(8*np.log(2))
                limits = (x0-3*sigma, x0+3*sigma)
            elif self.userot.get():
                x0 = hea[self.key_rotrv[0]]
                vsini = hea[self.key_rotvsini[0]]
                limits = (x0-vsini, x0+vsini)
            else:
                limits = (self.limitlow.get(), self.limitup.get())
            bis = pa.bisector(vrad, flux, limits=limits)
            bis_dict = {'rv_range' : bis['bisvel'], 'bisector' : bis['bisflux'], 'error' : bis['biserr']}
            self.update_text(f"Bisector:\n")
            self.update_text(f"Bis span = {np.round(bis['bispan'],4)} +/- {np.round(bis['e_bispan'],4)} km/s\n\n", timestamp=False)

            hea = self.insert_key(hea, self.key_ver[0], self.version, self.key_ver[1])   
            hea = self.insert_key(hea, self.pa_ver[0], self.pa_version, self.pa_ver[1])
            hea = self.insert_key(hea, self.key_lainp[0], os.path.basename(spec), self.key_lainp[1])                        
            hea = self.insert_key(hea, self.key_bislim[0], str(limits), self.key_bislim[1])
            hea = self.insert_key(hea, self.key_bispan[0], bis['bispan'], self.key_bispan[1])
            hea = self.insert_key(hea, self.key_biserr[0], bis['e_bispan'], self.key_biserr[1])
            self.save_fits(bis_dict,hea,spec,self.la_outdir.get(),self.bis_output.get())
            self.update_fits(spec,hea)
            
            if self.bisplot.get():
                self.xvalues = bis_dict['rv_range']
                self.yvalues = bis_dict['bisector']
                self.y_add = None
                self.y_res = None
                self.ymin = None
                self.ymax = None
                self.limits = (None,None)
                self.logscale = False
                self.do_plot = True


        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  END computing the bisector \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        #self.thread_finished()

        return

    def call_bisector(self):
        try:
            self.bisector()
        except:
            self.update_text("\n\n########################\nAborting the process.\nAre there any profiles with the same extension but different format?\n########################\n\n")
        self.thread_finished()        
        return

    def thread_bisector(self):
        # this thread call the "bisector" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_bisector = Thread(target=self.call_bisector, daemon=True)
            self.thread_running.set(1)
            self.th_bisector.start()
        return
         
    def fourier(self):
        # Compute Fourier transform
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  START computing the Fourier Transform \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        self.update_text("Parameters:\n", timestamp=False)
        self.update_text(f"Use Gaussian fitting: {self.usegauss.get()},\n", timestamp=False)
        self.update_text(f"Use rotational fitting: {self.userot.get()},\n", timestamp=False)
        self.update_text(f"Manual limits: ({self.limitlow.get()},{self.limitup.get()})  km/s\n", timestamp=False)
        self.update_text(f"Linear Limb Darkening: {self.ld.get()},\n", timestamp=False)

        spectra = self.find_file(self.la_spec.get(), self.la_indir.get())
        n_spectra = len(spectra)
        if not n_spectra:
            self.update_text(f"No input found. Aborted.\n")
            self.update_text("\n###########################\n", timestamp=False)
            self.update_text("  END computing the Fourier Transform \n", timestamp=False)
            self.update_text("###########################\n\n", timestamp=False)
            return
        self.update_text(f"Found {len(spectra)} profiles.\n")
        
        for n, spec in enumerate(spectra):
            if self.abort_value:
                self.update_text("Fourier Transform aborted.\n")
                self.abort_value = False
                break
            self.update_text(f"{n+1} of {n_spectra} profiles: working on profile {os.path.basename(spec)}.\n")
            with fits.open(spec) as hdu:
                data = hdu[1].data
                hea = hdu[0].header
            vrad = data.field(0)
            flux = data.field(1)
            try:
                errs = data.field(2)
            except KeyError:
                errs = False
            if self.usegauss.get():
                x0 = hea[self.key_gaussrv[0]]
                fwhm = hea[self.key_gaussfwhm[0]]
                sigma = fwhm/np.sqrt(8*np.log(2))
                limits = (x0-3*sigma, x0+3*sigma)
            elif self.userot.get():
                x0 = hea[self.key_rotrv[0]]
                vsini = hea[self.key_rotvsini[0]]
                limits = (x0-vsini, x0+vsini)
            else:
                limits = (self.limitlow.get(), self.limitup.get())
            fou = pa.fourier(vrad, flux, errs=errs, limits=limits, ld=self.ld.get())
            fou_dict = {'FFT_fr' : fou['FFT_fr'], 'FFT_pow' : fou['FFT_pow']}

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
            
            self.save_fits(fou_dict,hea,spec,self.la_outdir.get(),self.fou_output.get())
            self.update_fits(spec,hea)

            if self.fouplot.get():

                xmax = fou['zeros'][2]+0.5*fou['zeros'][2]
                xfind = np.searchsorted(fou_dict['FFT_fr'], xmax)
                yfind = np.searchsorted(fou_dict['FFT_fr'], 0.9*fou['zeros'][0], side='right')
                ymax = fou_dict['FFT_pow'][yfind]
                self.xvalues = fou_dict['FFT_fr'][:xfind]
                self.yvalues = fou_dict['FFT_pow'][:xfind]
                self.y_add = [np.ones(fou_dict['FFT_fr'][:xfind].shape)*fou['FFT_err']]
                self.y_res = None
                self.ymax = None
                self.ymin = None
                self.limits = (fou['zeros'][0],fou['zeros'][1],fou['zeros'][2])
                self.logscale = True
                self.do_plot = True

            self.update_text(f"Fourier Transform:\n")
            self.update_text(f"First zero vsini = {np.round(fou['vsini'][0],4)} +/- {np.round(fou['e_vsini'][0],4)} km/s\n", timestamp=False)
            self.update_text(f"Mean vsini = {np.round(fou['mean_vsini'],4)} +/- {np.round(fou['e_mean_vsini'],4)} km/s\n", timestamp=False)
            self.update_text(f"q2/q1 = {np.round(fou['ratio'],4)} +/- {np.round(fou['e_ratio'],4)} km/s\n\n", timestamp=False)
            self.update_text(f"RV = {np.round(fou['rv'],4)} km/s\n\n", timestamp=False)
            
        self.update_text("\n###########################\n", timestamp=False)
        self.update_text("  END computing the Fourier Transform \n", timestamp=False)
        self.update_text("###########################\n\n", timestamp=False)
        #self.thread_finished()

        return

    def call_fourier(self):
        try:
            self.fourier()
        except:
            self.update_text("\n\n########################\nAborting the process.\nAre there any profiles with the same extension but different format?\n########################\n\n")
        self.thread_finished()        
        return

    def thread_fourier(self):
        # this thread call the "fourier" function
        if self.thread_running.get():
            self.update_text("\n#######################\nAnother process is already running, wait for it to finish\n#######################\n", timestamp=False)
        else:
            self.th_fourier = Thread(target=self.call_fourier, daemon=True)
            self.thread_running.set(1)
            self.th_fourier.start()
        return

         
#####################################
if __name__ == '__main__':
   root = tk.Tk()
   main_app =  MainApp(root)

   root.mainloop()


