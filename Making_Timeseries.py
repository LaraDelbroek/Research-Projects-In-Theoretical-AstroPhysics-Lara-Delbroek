from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math
from scipy import interpolate
from scipy.optimize import curve_fit
import sys

# Original file by Nicolas Moens

bf = 'WR_2D_alpha_LTE_G4_O3_'
fn = bf + '0016.dat'

# Mean density
########################################
A = np.loadtxt(bf + '_rho', skiprows=1)
mean_r = np.transpose(A)[0]
mean_rho = np.transpose(A)[2]
f_rho = interpolate.interp1d(mean_r, mean_rho, kind='linear', fill_value="extrapolate")


snaps_arr = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 60, 80, 100]

fig, ax = plt.subplots(2, len(snaps_arr))
# fig.set_size_inches(10, 8)

for ii in range(len(snaps_arr)):
    # Retrieve data and coordinates from file
    fn = bf + str(snaps_arr[ii]).zfill(4) + '.dat'
    print(ii, ':', fn)
    ds = amrvac_reader.load_file(fn)
    ad = ds.load_all_data()
    x, y = ds.get_coordinate_arrays()
    t = ds.get_time()
    rho = ad['rho']
    v = ad['v1']

    relrho = rho.copy()

    for i in range(len(x)):
        for j in range(len(y)):
            mean_rho_local = f_rho(x[i])
            relrho[i, j] = rho[i, j] / mean_rho_local

    pcp0 = ax[0, ii].pcolor(y, x, relrho, norm=colors.LogNorm(vmin=0.1, vmax=10.0), cmap='RdYlBu')
    pcp1 = ax[1, ii].pcolor(y, x, v, norm=colors.TwoSlopeNorm(vmin=-0.25, vcenter=0., vmax=2.5), cmap='RdBu')

    ax[0, ii].set_aspect(aspect='equal')
    ax[1, ii].set_aspect(aspect='equal')

    ax[0, ii].set_title(str(t), rotation=45, y=1.0, pad=+15)
    ax[0, ii].set_yticks([])
    ax[1, ii].set_yticks([])
    ax[0, ii].set_xticks([])
    ax[1, ii].set_xticks([])

    plt.savefig('timeseries_rho_v.png')

ax[0, 0].set_title('$t/\\tau_{dyn} =$ \n 0 ', rotation=0, y=1.0, pad=+14)

ax[1, 0].set_yticks([1, 2, 3, 4, 5, 6])
ax[0, 0].set_yticks([1, 2, 3, 4, 5, 6])
ax[0, 0].set_ylabel("r $[R_{\\rm c}]$")
ax[1, 0].set_ylabel("r $[R_{\\rm c}]$")
ax[1, 0].set_xlabel("y $[R_{\\rm c}]$")

# cb0 = fig.colorbar(pcp0, ax=ax[0,:].ravel().tolist(),shrink=0.7)
# cb1 = fig.colorbar(pcp1, ax=ax[1,:].ravel().tolist(),shrink=0.7,ticks=[-0.25, 0, 2.5])

# cb0.set_label('relative density', rotation=270, fontsize=10, labelpad=20)
# cb1.set_label('radial velocity [1000 km $\\rm s^{-1}$]', rotation=270, fontsize=10, labelpad=20)

plt.tight_layout(h_pad=0.1, w_pad=0.0)

cb0 = fig.colorbar(pcp0, ax=ax[0, :].ravel().tolist(), shrink=0.7)
cb1 = fig.colorbar(pcp1, ax=ax[1, :].ravel().tolist(), shrink=0.7, ticks=[-0.25, 0, 2.5])

cb0.set_label('relative density', rotation=270, fontsize=10, labelpad=20)
cb1.set_label('radial velocity [1000 km $\\rm s^{-1}$]', rotation=270, fontsize=10, labelpad=20)

plt.savefig('timeseries_rot.png', bbox_inches='tight')
plt.show()
