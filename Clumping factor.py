import random
import statistics
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

fn = bf + '0016.dat'

ds = amrvac_reader.load_file(fn)  # Dit leest de data in de file
ad = ds.load_all_data()  # Dit voegt alle data toe in ad

x, y = ds.get_coordinate_arrays()  # Dit geeft afmetingen ofzo denk ik
# x gaat van 1 tot 6 in 1024 stappen. Dus dit lijkt hoogte te zijn.
# y gaat van -0.25 tot 0.25 in 128 stappen.
# Dus dit lijkt dan breedte te zijn, stappen van (ongeveer) 0.004

Height_tot = np.empty((len(Timestep_index), len(x)))
f_cl_tot = np.empty((len(Timestep_index), len(x)))

f_cl_tot_average = np.empty(len(x))

#######################################
# CLUMPING FACTOR
######################################"

for Timestep in range(0, len(Timestep_index)):
    ##################################
    # DATA RETRIEVAL
    ##################################

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

    ####################################
    # START CALCULATING
    ####################################

    # Formula= f_cl = average(rho^2) / (average(rho))^2

fn = bf + str(Timestep_index[0]).zfill(4) + '.dat'
ds = amrvac_reader.load_file(fn)
ad = ds.load_all_data()
x, y = ds.get_coordinate_arrays()
Combined_rho = ad['rho']
for index in range(1, len(Timestep_index)):  # MAKING ONE BIG BOX OUT OF THE RESULTS  STARTING FROM FILE 30
    fn = bf + str(Timestep_index[index]).zfill(4) + '.dat'
    ds = amrvac_reader.load_file(fn)
    ad = ds.load_all_data()
    x, y = ds.get_coordinate_arrays()
    rho = ad['rho']
    New_array = np.concatenate((Combined_rho, rho), axis=1)
    Combined_rho = New_array

rho_average = Combined_rho.sum(axis=1)/(128*len(Timestep_index))
rho_squared = np.square(Combined_rho)
rho_squared_average = rho_squared.sum(axis=1)/(128*len(Timestep_index))

f_cl =np.divide(rho_squared_average, np.square(rho_average))

fcl_flattened_per5 = []
for iii in range(0, len(f_cl) - 4, 5):
    fcl_flattened_per5.append(statistics.median([f_cl[iii], f_cl[iii + 1], f_cl[iii+ 2], f_cl[iii+ 3], f_cl[iii+ 4]]))

fcl_flattened_per10 = []
for iii in range(0, len(f_cl) - 9, 10):
    fcl_flattened_per10.append(
        statistics.median([f_cl[iii], f_cl[iii + 1], f_cl[iii + 2], f_cl[iii + 3], f_cl[iii + 4], f_cl[iii + 5],
        f_cl[iii + 6], f_cl[iii + 7], f_cl[iii + 8], f_cl[iii + 9]]))

fcl_flattened_per15 = []
for iii in range(0, len(f_cl) - 14, 15):
    fcl_flattened_per15.append(
        statistics.median([f_cl[iii], f_cl[iii + 1], f_cl[iii + 2], f_cl[iii + 3], f_cl[iii + 4], f_cl[iii + 5],
        f_cl[iii + 6], f_cl[iii + 7], f_cl[iii + 8], f_cl[iii + 9], f_cl[iii + 10], f_cl[iii + 11], f_cl[iii + 12],
         f_cl[iii + 13], f_cl[iii + 14]]))

fcl_flattened_per20 = []
for iii in range(0, len(f_cl) - 19, 20):
    fcl_flattened_per20.append(
        statistics.median([f_cl[iii], f_cl[iii + 1], f_cl[iii + 2], f_cl[iii + 3], f_cl[iii + 4], f_cl[iii + 5],
        f_cl[iii + 6], f_cl[iii + 7], f_cl[iii + 8], f_cl[iii + 9], f_cl[iii + 10], f_cl[iii + 11], f_cl[iii + 12],
         f_cl[iii + 13], f_cl[iii + 14], f_cl[iii + 15], f_cl[iii + 16], f_cl[iii + 17], f_cl[iii + 18],
         f_cl[iii + 19]]))

fcl_flattened_per50 = []
for iii in range(0, len(f_cl) - 49, 50):
    fcl_flattened_per50.append(
        statistics.median([f_cl[iii], f_cl[iii + 1], f_cl[iii + 2], f_cl[iii + 3], f_cl[iii + 4], f_cl[iii + 5],
        f_cl[iii + 6], f_cl[iii + 7], f_cl[iii + 8], f_cl[iii + 9], f_cl[iii + 10], f_cl[iii + 11], f_cl[iii + 12],
         f_cl[iii + 13], f_cl[iii + 14], f_cl[iii + 15], f_cl[iii + 16], f_cl[iii + 17], f_cl[iii + 18],
         f_cl[iii + 19], f_cl[iii + 20], f_cl[iii + 21], f_cl[iii + 22], f_cl[iii + 23], f_cl[iii + 24],
         f_cl[iii + 25], f_cl[iii + 26], f_cl[iii + 27], f_cl[iii + 28], f_cl[iii + 29], f_cl[iii + 30],
         f_cl[iii + 31], f_cl[iii + 32], f_cl[iii + 33], f_cl[iii + 34], f_cl[iii + 35], f_cl[iii + 36],
         f_cl[iii + 37], f_cl[iii + 38], f_cl[iii + 39], f_cl[iii + 40], f_cl[iii + 41], f_cl[iii + 42],
         f_cl[iii + 43], f_cl[iii + 44], f_cl[iii + 45], f_cl[iii + 46], f_cl[iii + 47], f_cl[iii + 48],
         f_cl[iii + 49]]))



Height = []
for height in range(0, len(x)):  # Opstellen hoogte voor plotten
    Height.append(x[height])

Height_flattened_per5 =  []
for height in range(0, len(x) - 4, 5):
    Height_flattened_per5.append(x[height])

Height_flattened_per10 =  []
for height in range(0, len(x) - 9, 10):
    Height_flattened_per10.append(x[height])

Height_flattened_per15 =  []
for height in range(0, len(x) - 14, 15):
    Height_flattened_per15.append(x[height])

Height_flattened_per20 =  []
for height in range(0, len(x) - 19, 20):
    Height_flattened_per20.append(x[height])

Height_flattened_per50 =  []
for height in range(0, len(x) - 49, 50):
    Height_flattened_per50.append(x[height])

bf = 'WR_2D_alpha_LTE_G4_O3_'
Timestep_index = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
                  55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                  80, 71, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
fn = bf + '0016.dat'

ds = amrvac_reader.load_file(fn)  # Dit leest de data in de file
ad = ds.load_all_data()  # Dit voegt alle data toe in ad

x, y = ds.get_coordinate_arrays()  # Dit geeft afmetingen ofzo denk ik
# x gaat van 1 tot 6 in 1024 stappen. Dus dit lijkt hoogte te zijn.
# y gaat van -0.25 tot 0.25 in 128 stappen.
# Dus dit lijkt dan breedte te zijn, stappen van (ongeveer) 0.004

Height_tot = np.empty((len(Timestep_index), len(x)))
f_cl_tot = np.empty((len(Timestep_index), len(x)))

f_cl_tot_average = np.empty(len(x))

#######################################
# CLUMPING FACTOR
######################################"

for Timestep in range(0, len(Timestep_index)):
    ##################################
    # DATA RETRIEVAL
    ##################################

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

    ####################################
    # START CALCULATING
    ####################################

    # Formula= f_cl = average(rho^2) / (average(rho))^2

fn = bf + str(Timestep_index[0]).zfill(4) + '.dat'
ds = amrvac_reader.load_file(fn)
ad = ds.load_all_data()
x, y = ds.get_coordinate_arrays()
Combined_rho = ad['rho']
for index in range(1, len(Timestep_index)):  # MAKING ONE BIG BOX OUT OF THE RESULTS  STARTING FROM FILE 30
    fn = bf + str(Timestep_index[index]).zfill(4) + '.dat'
    ds = amrvac_reader.load_file(fn)
    ad = ds.load_all_data()
    x, y = ds.get_coordinate_arrays()
    rho = ad['rho']
    New_array = np.concatenate((Combined_rho, rho), axis=1)
    Combined_rho = New_array

rho_average = Combined_rho.sum(axis=1)/(128*len(Timestep_index))
rho_squared = np.square(Combined_rho)
rho_squared_average = rho_squared.sum(axis=1)/(128*len(Timestep_index))

f_cl_rot =np.divide(rho_squared_average, np.square(rho_average))

fcl_flattened_per5_rot = []
for iii in range(0, len(f_cl) - 4, 5):
    fcl_flattened_per5_rot.append(statistics.median([f_cl_rot[iii], f_cl_rot[iii + 1], f_cl_rot[iii+ 2], f_cl_rot[iii+ 3], f_cl_rot[iii+ 4]]))

fcl_flattened_per10_rot = []
for iii in range(0, len(f_cl) - 9, 10):
    fcl_flattened_per10_rot.append(
        statistics.median([f_cl_rot[iii], f_cl_rot[iii + 1], f_cl_rot[iii + 2], f_cl_rot[iii + 3], f_cl_rot[iii + 4], f_cl_rot[iii + 5],
        f_cl_rot[iii + 6], f_cl_rot[iii + 7], f_cl_rot[iii + 8], f_cl_rot[iii + 9]]))

fcl_flattened_per15_rot = []
for iii in range(0, len(f_cl) - 14, 15):
    fcl_flattened_per15_rot.append(
        statistics.median([f_cl_rot[iii], f_cl_rot[iii + 1], f_cl_rot[iii + 2], f_cl_rot[iii + 3], f_cl_rot[iii + 4], f_cl_rot[iii + 5],
        f_cl_rot[iii + 6], f_cl_rot[iii + 7], f_cl_rot[iii + 8], f_cl_rot[iii + 9], f_cl_rot[iii + 10], f_cl_rot[iii + 11], f_cl_rot[iii + 12],
         f_cl_rot[iii + 13], f_cl_rot[iii + 14]]))

fcl_flattened_per20_rot = []
for iii in range(0, len(f_cl) - 19, 20):
    fcl_flattened_per20_rot.append(
        statistics.median([f_cl_rot[iii], f_cl_rot[iii + 1], f_cl_rot[iii + 2], f_cl_rot[iii + 3], f_cl_rot[iii + 4], f_cl_rot[iii + 5],
        f_cl_rot[iii + 6], f_cl_rot[iii + 7], f_cl_rot[iii + 8], f_cl_rot[iii + 9], f_cl_rot[iii + 10], f_cl_rot[iii + 11], f_cl_rot[iii + 12],
         f_cl_rot[iii + 13], f_cl_rot[iii + 14], f_cl_rot[iii + 15], f_cl_rot[iii + 16], f_cl_rot[iii + 17], f_cl_rot[iii + 18],
         f_cl_rot[iii + 19]]))

fcl_flattened_per50_rot = []
for iii in range(0, len(f_cl) - 49, 50):
    fcl_flattened_per50_rot.append(
        statistics.median([f_cl_rot[iii], f_cl_rot[iii + 1], f_cl_rot[iii + 2], f_cl_rot[iii + 3], f_cl_rot[iii + 4], f_cl_rot[iii + 5],
        f_cl_rot[iii + 6], f_cl_rot[iii + 7], f_cl_rot[iii + 8], f_cl_rot[iii + 9], f_cl_rot[iii + 10], f_cl_rot[iii + 11], f_cl_rot[iii + 12],
         f_cl_rot[iii + 13], f_cl_rot[iii + 14], f_cl_rot[iii + 15], f_cl_rot[iii + 16], f_cl_rot[iii + 17], f_cl_rot[iii + 18],
         f_cl_rot[iii + 19], f_cl_rot[iii + 20], f_cl_rot[iii + 21], f_cl_rot[iii + 22], f_cl_rot[iii + 23], f_cl_rot[iii + 24],
         f_cl_rot[iii + 25], f_cl_rot[iii + 26], f_cl_rot[iii + 27], f_cl_rot[iii + 28], f_cl_rot[iii + 29], f_cl_rot[iii + 30],
         f_cl_rot[iii + 31], f_cl_rot[iii + 32], f_cl_rot[iii + 33], f_cl_rot[iii + 34], f_cl_rot[iii + 35], f_cl_rot[iii + 36],
         f_cl_rot[iii + 37], f_cl_rot[iii + 38], f_cl_rot[iii + 39], f_cl_rot[iii + 40], f_cl_rot[iii + 41], f_cl_rot[iii + 42],
         f_cl_rot[iii + 43], f_cl_rot[iii + 44], f_cl_rot[iii + 45], f_cl_rot[iii + 46], f_cl_rot[iii + 47], f_cl_rot[iii + 48],
         f_cl_rot[iii + 49]]))



Height_rot = []
for height in range(0, len(x)):  # Opstellen hoogte voor plotten
    Height_rot.append(x[height])

Height_flattened_per5_rot =  []
for height in range(0, len(x) - 4, 5):
    Height_flattened_per5_rot.append(x[height])

Height_flattened_per10_rot =  []
for height in range(0, len(x) - 9, 10):
    Height_flattened_per10_rot.append(x[height])

Height_flattened_per15_rot =  []
for height in range(0, len(x) - 14, 15):
    Height_flattened_per15_rot.append(x[height])

Height_flattened_per20_rot =  []
for height in range(0, len(x) - 19, 20):
    Height_flattened_per20_rot.append(x[height])

Height_flattened_per50_rot =  []
for height in range(0, len(x) - 49, 50):
    Height_flattened_per50_rot.append(x[height])

plt.rcParams['font.size'] = 12
plt.plot(Height, f_cl, label ="No rotation")
plt.plot(Height_rot, f_cl_rot, label ="Rotation")  # PLOTTING AVERAGED RESULTS
plt.xlabel('Height (Rc)')
plt.ylabel('clumping factor')
plt.title('Clumping factor')
plt.legend(loc="upper right")
plt.savefig('Clumping11', bbox_inches='tight')
plt.figure()
plt.plot(Height_flattened_per5, fcl_flattened_per5, label ="No rotation")              # PLOTTING AVERAGED RESULTS
plt.plot(Height_flattened_per5_rot, fcl_flattened_per5_rot, label ="Rotation")
plt.xlabel('Height (Rc)')
plt.ylabel('clumping factor')
plt.title('Clumping factor flattened per 5 cells')
plt.legend(loc="upper right")
plt.savefig('Clumping22', bbox_inches='tight')
plt.figure()
plt.plot(Height_flattened_per10, fcl_flattened_per10, label ="No rotation")
plt.plot(Height_flattened_per10_rot, fcl_flattened_per10_rot, label ="Rotation")  # PLOTTING AVERAGED RESULTS
plt.xlabel('Height (Rc)')
plt.ylabel('clumping factor')
plt.title('Clumping factor flattened per 10 cells')
plt.legend(loc="upper right")
plt.savefig('Clumping33', bbox_inches='tight')
plt.figure()
plt.plot(Height_flattened_per15, fcl_flattened_per15, label ="No rotation")
plt.plot(Height_flattened_per15_rot, fcl_flattened_per15_rot, label ="Rotation")  # PLOTTING AVERAGED RESULTS
plt.xlabel('Height (Rc)')
plt.ylabel('clumping factor')
plt.title('Clumping factor flattened per 15 cells')
plt.legend(loc="upper right")
plt.savefig('Clumping44', bbox_inches='tight')
plt.figure()
plt.plot(Height_flattened_per20, fcl_flattened_per20, label ="No rotation")
plt.plot(Height_flattened_per20_rot, fcl_flattened_per20_rot, label ="Rotation") # PLOTTING AVERAGED RESULTS
plt.xlabel('Height (Rc)')
plt.ylabel('clumping factor')
plt.title('Clumping factor flattened per 20 cells')
plt.legend(loc="upper right")
plt.savefig('Clumping55', bbox_inches='tight')
plt.figure()
plt.plot(Height_flattened_per50, fcl_flattened_per50, label ="No rotation")
plt.plot(Height_flattened_per50_rot, fcl_flattened_per50_rot, label ="Rotation") # PLOTTING AVERAGED RESULTS
plt.xlabel('Height (Rc)')
plt.ylabel('clumping factor')
plt.title('Clumping factor flattened per 50 cells')
plt.legend(loc="upper right")
plt.savefig('Clumping66', bbox_inches='tight')
plt.show()

