##############################################################################################
# INITIALIZE PACKAGES
##############################################################################################

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import FK5
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import healpy as hp
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy import units
from scipy.interpolate import RectBivariateSpline
from scipy import interpolate
import copy
import os
import sys
from tqdm import tqdm
import glob
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams["font.size"] = "18"

import time
import pandas as pd
start = time.time()


##############################################################################################
# SPECIFY THE PATH WHERE THE ASSOCIATED FILES ARE PRESENT
##############################################################################################

PATH="/home/abhijeet/Saras-lcrt/modified_files"
sys.path.insert(1, PATH)

from observe import Observations # This calls all the classes and functions in the files specified in PATH

##############################################################################################
# LCRT and SARAS
##############################################################################################

lcrt = False
lcrt_with_cos_beam = True
saras = False
def get_freq_from_lcrt(filename):  # This function simply reads the number of the filename e.g. "ticra-export_5MHz.dat" ---> [5.0]
    _temp = os.path.basename(filename).replace('ticra-export_','').replace('MHz.dat','')
    return float(_temp)

if lcrt == True:
    bm = False                                 # We don't use beam files of saras
    gm = False           
    SITE_LATITUDE   = 0.0
    SITE_LONGITUDE  = 180.0
    ELEVATION       = 0.0
    
    fmin = 11/1e3
    fmax = 47/1e3

    phi_res   = 5.0
    theta_res = 0.2

    theta_array = np.arange(-90, 90 + theta_res, theta_res)
    theta_array = np.round(theta_array, 2)    # For multiples of 0.6, the values were 0.600001. Don't know the reason
    phi_array = np.arange(0, 360 + phi_res, phi_res)
    t1 = []
    for i in range(len(phi_array)):
        x = phi_array[i]
        t = list(np.tile(x,901))
        t1.extend(t)

    tt = np.array(t1)
    phi_col =  tt

    lcrt_path = '/home/abhijeet/Saras-lcrt/lcrt_reflector_patterns'
    lcrt_list = []
    lcrt_list = sorted(glob.glob(os.path.join(lcrt_path, '*.dat')))
    freq_array  = []
    for ii, file_add in enumerate(lcrt_list):
        freq_array.append(get_freq_from_lcrt(file_add)) # This is simply the array of numbers of the files names

        
    isort   = np.argsort(np.array(freq_array))
    freq_lcrt    = np.array(freq_array)[isort]
    freq_lcrt    = freq_lcrt[freq_lcrt > 10]
    dt = TimeDelta(np.linspace(0.,655.2*3600, 2), format='sec')

elif lcrt_with_cos_beam == True:
    bm = False          # We don't use beam files of saras
    gm = False       
    SITE_LATITUDE   = 0.0
    SITE_LONGITUDE  = 180
    ELEVATION       = 0.0
    
    fmin = 11/1e3
    fmax = 47/1e3

    phi_res   = 5.0
    theta_res = 0.2

    theta_array = np.arange(-90, 90 + theta_res, theta_res)
    theta_array = np.round(theta_array, 2)    # For multiples of 0.6, the values were 0.600001. Don't know the reason
    phi_array = np.arange(0, 360 + phi_res, phi_res)
    t1 = []
    for i in range(len(phi_array)):
        x = phi_array[i]
        t = list(np.tile(x,901))
        t1.extend(t)

    tt = np.array(t1)
    phi_col =  tt

    lcrt_path = '/home/abhijeet/Saras-lcrt/lcrt_reflector_patterns'
    lcrt_list = []
    lcrt_list = sorted(glob.glob(os.path.join(lcrt_path, '*.dat')))
    freq_array  = []
    for ii, file_add in enumerate(lcrt_list):
        freq_array.append(get_freq_from_lcrt(file_add)) # This is simply the array of numbers of the files names

    isort   = np.argsort(np.array(freq_array))
    freq_lcrt    = np.array(freq_array)[isort]
    freq_lcrt    = freq_lcrt[freq_lcrt > 10]
    dt = TimeDelta(np.linspace(0.,655.2*3600, 2), format='sec')   

else:
    bm = True
    gm = True
    SITE_LATITUDE       = 32.81302 #/* +32.77944444   deg for Hanle */ 
    SITE_LONGITUDE      = 78.87130 #/* 78.96416667 deg for Hanle */
    ELEVATION = 4500

    #fmin = 56/1e3
    #fmax = 109/1e3
    fmin = 55/1e3
    fmax = 110/1e3
    phi_res   = 1
    theta_res = 1

    theta_array = np.arange(-90, 90 + theta_res, theta_res)       
    phi_array   = np.arange(0, 360, phi_res)
    dt = TimeDelta(np.linspace(0.,24.*3600, 10), format='sec')

freq_width = 2 # For model without data files

obstimes = Time('2019-4-12 23:00:00') + dt

##############################################################################################
# SELECT FILES OR THE FUNCTION
##############################################################################################

beam_path = '/home/abhijeet/Saras-codes-main/elec_small_ant_with_2_12U_bus_cone_scale_0_925_ref_scale_0_9999'
lcrt_path = '/home/abhijeet/Saras-lcrt/lcrt_reflector_patterns'

beam_saras = []
rt_file = []
beam_lcrt = []

if saras == True:
    bm == True
    gm == True
    
if lcrt == True:
    beam_lcrt = sorted(glob.glob(os.path.join(lcrt_path, '*.dat')))
    
if bm == True:
    beam_saras = sorted(glob.glob(os.path.join(beam_path, '*mhz.txt')))
else:
    beam_saras = []

if gm == True:
    rt_file = os.path.join(beam_path, 'gamma_linear.txt')
elif lcrt == True:
    rt_file =  os.path.join(lcrt_path, 's11_new_LPDA_feed.txt')
else:
    rt_file = []


output_name = "fiducial_monopole"

##############################################################################################
# RUN COMMAND
##############################################################################################

param_sky = {'fmin': fmin, 'fmax': fmax}

if saras == True:
    param_obs = {'file_saras': beam_saras,\
                 'gamma_file': rt_file,\
                 'phi_array': phi_array,\
                 'theta_array': theta_array,\
                 'fmin':fmin,\
                 'fmax':fmax,\
                 'SITE_LATITUDE':SITE_LATITUDE,\
                 'SITE_LONGITUDE':SITE_LONGITUDE,\
                 'ELEVATION':ELEVATION,\
                 'obstimes':obstimes,\
                 'include_gamma':True,\
                 'output':True,\
                 'file_name': output_name,\
                 'path_output': '/home/abhijeet/Saras-lcrt/output'}
    
elif lcrt == True:                      # Here we are including the frequency file in param_obs
    param_obs = {'file_lcrt': beam_lcrt,\
                 'gamma_file': rt_file,\
                 'cos_beam': False,\
                 'freq': freq_lcrt,\
                 'phi_col': phi_col,\
                 'phi_array': phi_array,\
                 'theta_array': theta_array,\
                 'fmin':fmin,\
                 'fmax':fmax,\
                 'SITE_LATITUDE':SITE_LATITUDE,\
                 'SITE_LONGITUDE':SITE_LONGITUDE,\
                 'ELEVATION':ELEVATION,\
                 'obstimes':obstimes,\
                 'include_gamma':False,\
                 'output':True,\
                 'file_name': output_name,\
                 'path_output': '/home/abhijeet/Saras-lcrt/output'}
    
elif lcrt_with_cos_beam == True:
    param_obs = {'file_lcrt': beam_lcrt,\
                 'cos_beam': True,\
                 'freq': freq_lcrt,\
                 'phi_col': phi_col,\
                 'phi_array': phi_array,\
                 'theta_array': theta_array,\
                 'fmin':fmin,\
                 'fmax':fmax,\
                 'SITE_LATITUDE':SITE_LATITUDE,\
                 'SITE_LONGITUDE':SITE_LONGITUDE,\
                 'ELEVATION':ELEVATION,\
                 'obstimes':obstimes,\
                 'include_gamma':False,\
                 'output':True,\
                 'file_name': output_name,\
                 'path_output': '/home/abhijeet/Saras-lcrt/output'}
    
else:
    param_obs = {'phi_array': phi_array,\
                 'theta_array': theta_array,\
                 'fmin':fmin,\
                 'fmax':fmax,\
                 'freq_width':freq_width,\
                 'SITE_LATITUDE':SITE_LATITUDE,\
                 'SITE_LONGITUDE':SITE_LONGITUDE,\
                 'ELEVATION':ELEVATION,\
                 'obstimes':obstimes,\
                 'include_gamma':True,\
                 'output':True,\
                 'file_name': output_name,\
                 'path_output': '/home/abhijeet/Saras-lcrt/output'}



TA = Observations(**param_obs)
TA.make_observations()

end = time.time()
print("Time taken:", end - start, "seconds")