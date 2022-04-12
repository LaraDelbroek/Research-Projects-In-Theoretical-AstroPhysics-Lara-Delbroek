import random

from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math

from scipy import interpolate
from scipy.optimize import curve_fit
import sys

############################################
# DEFINING NEEDED THINGS SO THE SCRIPT WORKS
#############################################

bf = 'WR_Isotropic_Calc1alpha_LTE_'
fn = bf + '0016.dat'

Timestep_index = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
                  55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                  80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                  104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 115, 116, 117, 118, 119, 120, 121, 122, 123,
                  124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
                  144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163,
                  164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                  184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]

# Nummer 114 weggehaald omdat er een numerieke fout was.

ds = amrvac_reader.load_file(fn)  # Dit leest de data in de file
ad = ds.load_all_data()  # Dit voegt alle data toe in ad

x, y = ds.get_coordinate_arrays()  # Dit geeft afmetingen ofzo denk ik
# x gaat van 1 tot 6 in 1024 stappen. Dus dit lijkt hoogte te zijn.
# y gaat van -0.25 tot 0.25 in 128 stappen.
# Dus dit lijkt dan breedte te zijn, stappen van (ongeveer) 0.004

v_r_disp_tot = np.empty((len(Timestep_index), len(x)))
v_t_disp_tot = np.empty((len(Timestep_index), len(x)))
v_t_RMS_tot = np.empty((len(Timestep_index), len(x)))
Height_tot = np.empty((len(Timestep_index), len(x)))
f_cl_tot = np.empty((len(Timestep_index), len(x)))

v_r_disp_tot_average = np.empty(len(x))
v_t_disp_tot_average = np.empty(len(x))
f_cl_tot_average = np.empty(len(x))

###############################
# v_disp
###############################

for Timestep in range(0, len(Timestep_index)):

    ##################################
    # DATA RETRIEVAL
    ##################################

    # Formula: v_r_disp = sqrt[(average (v_r^2) - (average v_r)^2]
    # Formula: v_t_disp = sqrt[average (v_t^2) - (average v_t)^2]
    # v_t_disp should be equal to sqrt[average (v_t^2)] due to no net displacement expected in t direction

    fn = bf + str(Timestep_index[Timestep]).zfill(4) + '.dat'
    print("File 0", ':', fn)
    ds = amrvac_reader.load_file(fn)  # Dit leest de data in de file
    ad = ds.load_all_data()  # Dit voegt alle data toe in ad

    x, y = ds.get_coordinate_arrays()  # Dit geeft afmetingen ofzo denk ik
    # x gaat van 1 tot 6 in 1024 stappen. Dus dit lijkt hoogte te zijn.
    # y gaat van -0.25 tot 0.25 in 128 stappen.
    # Dus dit lijkt dan breedte te zijn, stappen van (ongeveer) 0.004

    t = ds.get_time()  # Dit geeft tijdstappen weer, specifiek de tijdstap
    # (of eenheid) waar de calculatie in zit

    rho = ad['rho']  # Dit geeft de rho data
    # Er zijn 1024 rijen en 128 kolommen, dus elke rij is de data van 1 r waarde

    v_r = ad['v1']  # Dit geeft de radiale snelheid
    # Analoge verdeling als bij rho
    v_t = ad['v2']  # Transversale snelheid

    ##################################
    # START CALCULATING
    ###################################

    fn = bf + str(Timestep_index[0]).zfill(4) + '.dat'
    ds = amrvac_reader.load_file(fn)
    ad = ds.load_all_data()
    x, y = ds.get_coordinate_arrays()
    Combined_v_r = ad['v1']
    Combined_v_t = ad['v2']

for index in range(1, len(Timestep_index)):  # MAKING ONE BIG BOX OUT OF THE RESULTS  STARTING FROM FILE 30
    fn = bf + str(Timestep_index[index]).zfill(4) + '.dat'
    ds = amrvac_reader.load_file(fn)
    ad = ds.load_all_data()
    x, y = ds.get_coordinate_arrays()
    vr = ad['v1']
    vt = ad['v2']
    New_array1 = np.concatenate((Combined_v_r, vr), axis=1)
    New_array2 = np.concatenate((Combined_v_t, vt), axis=1)
    Combined_v_r = New_array1
    Combined_v_t = New_array2

v_r_sum = Combined_v_r.sum(axis=1)  # Som van alle cellen per hoogte
v_r_squared = np.square(Combined_v_r)  # Squared
v_r_sum_squared = v_r_squared.sum(axis=1)  # Som van alle squared cellen per hoogte
v_t_sum = Combined_v_t.sum(axis=1)  # Som van alle cellen per hoogte
v_t_squared = np.square(Combined_v_t)  # Squared
v_t_sum_squared = v_t_squared.sum(axis=1)  # Som van alle squared cellen per hoogte

v_r_average = v_r_sum / (128 * len(Timestep_index))
v_r_squared_average = v_r_sum_squared / (128 * len(Timestep_index))
v_t_average = v_t_sum / (128 * len(Timestep_index))
v_t_squared_average = v_t_sum_squared / (128 * len(Timestep_index))

Height = []
for height in range(0, len(x)):  # Opstellen hoogte voor plotten
    Height.append(x[height])

    v_r_disp = np.sqrt(v_r_squared_average - np.square(v_r_average))  # USED FORMULAS
    v_t_disp = np.sqrt(v_t_squared_average - np.square(v_t_average))
    v_t_RMS = np.sqrt(v_t_squared_average)

#####################
# FINAL PLOTTING #
#####################
plt.rcParams['font.size'] = 16
plt.plot(Height, v_r_disp)  # AVERAGED RESULT
plt.xlabel('Height (Rc)')
plt.ylabel('dispersion(1000 km/s)')
plt.savefig('vr_disp_isotropic.png')
plt.show()

plt.rcParams['font.size'] = 18
plt.plot(Height, v_t_disp)  # AVERAGED RESULT
plt.xlabel('Height (Rc)')
plt.ylabel('dispersion(1000 km/s)')
plt.savefig('vt_disp_isotropic.png')
plt.show()

plt.rcParams['font.size'] = 18
plt.plot(Height, v_t_RMS)
plt.xlabel('Height (Rc)')
plt.ylabel('dispersion(1000 km/s)')
plt.savefig('vt_RMS_isotropic.png')
plt.show()
