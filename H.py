import random

from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math
from scipy.optimize import curve_fit


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
                  80, 71, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                  104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
                  124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
                  144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163,
                  164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                  184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]

ds = amrvac_reader.load_file(fn)  # Dit leest de data in de file
ad = ds.load_all_data()  # Dit voegt alle data toe in ad

x, y = ds.get_coordinate_arrays()  # Dit geeft afmetingen ofzo denk ik
# x gaat van 1 tot 6 in 1024 stappen. Dus dit lijkt hoogte te zijn.
# y gaat van -0.25 tot 0.25 in 128 stappen.
# Dus dit lijkt dan breedte te zijn, stappen van (ongeveer) 0.004

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

############################
# DETERMINING H
############################
H_array = []
delta_y_array = []
t0_array = []
print(np.mean(Combined_rho, axis=1))
for ii in range(len(x)):
    Combined_rho[ii,:] = Combined_rho[ii,:]/np.mean(Combined_rho[ii,:])
# Eerste massa per delta_y bepalen

for delta_y in np.logspace(2, 3.2, 20):   # is gaan van 10 tot 10000 in 100 logspace verdeelde stappen
    delta_y = int(delta_y)
    delta_y_array.append(delta_y)
    print(delta_y)
    Mass_distribution = []
    for Teller in range(0,20000):     # aantal keer dat m bepaald wordt op delta y
        Mass = 0
        y_coord = (Teller*13)   # Coordinaten bepalen  % modulo operator = delen en er dan int van maken
        x_coord = (Teller*7)%(307)
        for Stappen in range(0,delta_y):    # Massa bepalen
            Mass = Mass + Combined_rho[x_coord, (y_coord + Stappen)%(128*len(Timestep_index))]     # Hier nog toevoegen wat de breedte per cel is
        Mass_distribution.append(Mass/delta_y)


    ###############################
    # Fitting data to curve
    ##############################

    def fit_function(m, C, n, t0):
        return (C * (m ** n) * np.exp(-m * t0))


    bins = np.linspace(0, 4, 50)
    Massa, bins1 = np.histogram(Mass_distribution, bins=bins)
    bincenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])

    popt, pcov = curve_fit(fit_function, xdata=bincenters, ydata=Massa)
    print(popt)
    print(popt[2])
    xspace = np.linspace(0, 4, 100000)

    # t0 - 1 = delta_y/H
    H = delta_y/(popt[2] - 1)
    H_array.append(H)
    t0_array.append(popt[2])
    print(H)

    #plt.bar(bincenters, Massa, width=bins[1] - bins[0], color='navy', label=r'Histogram Mass distribution')
    #plt.plot(xspace, fit_function(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')

    # plt.xlim(0, 4)
    # plt.xlabel(r'Mass distribution')
    # plt.ylabel(r'Probability density function')
    # plt.title('Porositeitsfactor H =  %i' %H)
    # plt.legend(loc='best')
    # plt.show()
    # plt.clf()

def fit_function(delta_y, H):
    t0 = delta_y/H + 1
    return (t0)

popt, pcov = curve_fit(fit_function, xdata=delta_y_array, ydata=t0_array)
print(popt)
print('H waarde = ', popt[0])


plt.plot(delta_y_array, H_array)
plt.figure()
plt.plot(delta_y_array, t0_array)
plt.plot(delta_y_array, fit_function(delta_y_array, popt[0]))

plt.show()

