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
Timestep_index = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
                  55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                  80, 71, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                  104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
                  124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
                  144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163,
                  164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                  184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]


#####################################################
# DENSITY AUTOCORRELATION LENGTH
###############################################
# Formula to use:  fc(delta_y) = sum(t)sum(y) (rho(y) - average(rho))*(rho(y + delta_y) - average(rho))
# Delta_y is the thing that has to be varied here
# This for all r under 2.5 (index 307) and for all r above 2.5

fn = bf + str(Timestep_index[0]).zfill(4) + '.dat'
ds = amrvac_reader.load_file(fn)
ad = ds.load_all_data()
x, y = ds.get_coordinate_arrays()
Combined_rho = ad['rho']

for index in range (1, len(Timestep_index)):     # MAKING ONE BIG BOX OUT OF THE RESULTS  STARTING FROM FILE 30
    fn = bf + str(Timestep_index[index]).zfill(4) + '.dat'
    ds = amrvac_reader.load_file(fn)
    ad = ds.load_all_data()
    x, y = ds.get_coordinate_arrays()
    rho = ad['rho']
    New_array = np.concatenate((Combined_rho, rho), axis=1)
    Combined_rho = New_array


fc_delta_y = []
x_as = []
for Delta_y in range (-20, 20, 1):    # FOR A NUMBER OF DELTA_Y    # WEET NIET OF DIT OKE IS, VANAF WELKE WAARDEN DIT
    print(Delta_y)                                                  # EEN CORRECT BEELD GEEFT. DUURT LANG
    fc_altitude = 0
    x_as.append(Delta_y)
    for hoogte in range (100, 300,100): # UNDER 2.5: under 307
        fc_row = 0
        for ij in range (0 + abs(Delta_y) + 1, 128*len(Timestep_index) - abs(Delta_y) - 1): # SUMMATION OVER Y
            rho_average = Combined_rho.sum(axis=1) / (128 * len(Timestep_index))
            fc_lonely = ((Combined_rho[hoogte, ij] - rho_average[hoogte])*     # Voor 1 ij
                         (Combined_rho[hoogte, ij + Delta_y] - rho_average[hoogte]))
            fc_row = fc_row + fc_lonely       # Dus alle resultaten per rij worden opgeteld
        fc_altitude = fc_altitude + fc_row     # Per hoogte wordt het resultaat van een rij bijgeteld
    fc_delta_y.append(fc_altitude)

# FULL WIDTH HALF MAXIMUM BEPALEN
# Piek op delta_y = 0
Center = x_as.index(0)
Average_Height = 100000
Jump = 1

while Average_Height > (fc_delta_y[Center]/2):
    Average_Height = (fc_delta_y[Center - Jump] + fc_delta_y[Center + Jump])/2
    Jump = Jump + 1
    print(Jump)

Jump = Jump - 1
print(Jump)

FWHM = x_as[Center + Jump] - x_as[Center - Jump]

print('FWHM of fc curve = ',FWHM, 'cells')
print(Center)
print(Jump)

plt.plot(x_as, fc_delta_y)
plt.axvline(x = x_as[Center + Jump])
plt.axvline(x = x_as[Center - Jump])
plt.savefig('fc_isotropic.png')
plt.show()

np.savetxt("x_as_under_2_5.txt", x_as, delimiter=",")
np.savetxt("fc_delta_y_under_2_2.txt", fc_delta_y, delimiter=",")














