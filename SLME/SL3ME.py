import numpy as np
import csv
import copy
import random
import matplotlib.pyplot as plt
import pandas
import sklearn
from sklearn import metrics


material_direct_allowed_gap = 2.304
material_indirect_gap = 2.294




def calculate_SLME(material_eV_for_absorbance_data, material_absorbance_data,
                   material_direct_allowed_gap, material_indirect_gap,
                   thickness=50E-6, T=293.15, absorbance_in_inverse_centimeters=False,
                   cut_off_absorbance_below_direct_allowed_gap=True):
    ##### IMPORTANT NOTES:
    ## Material calculated absorbance is assumed to be in m^-1, not cm^-1!  (Most sources will provide absorbance in cm^-1, so be careful.)
    ## The default is to remove absorbance below the direct allowed gap. This is for dealing with broadening applied in DFT absorbance calculations. Probably not desired for experimental data.
    ## We can calculate at different temperatures if we want to, but 25 C / 293.15 K is the standard temperature assumed if not specified

    ## If absorbance is in cm^-1, multiply values by 100 to match units assumed in code
    if absorbance_in_inverse_centimeters:
        material_absorbance_data = material_absorbance_data * 100

    from numpy import sys, loadtxt, pi, exp, zeros, linspace
    # For flat plate solar panels, we want the "Global Tilt" spectra, this file is assumed to be in the directory
    try:
        solar_spectra_wavelength, solar_spectra_irradiance = loadtxt("am1.5G.dat", usecols=[0, 1], unpack=True,
                                                                     skiprows=2)
    except OSError:
        print('Could not locate am1.5G.dat datafile. Please place this file in the local directory.')
        sys.exit()
    solar_spectra_wavelength_meters = solar_spectra_wavelength * 10 ** -9

    c = 299792458  # speed of light, m/s
    h = 6.62607004081E-34  # Planck's constant J*s (W)
    h_eV = 4.135667516E-15  # Planck's constant eV*s
    k = 1.3806485279E-23  # Boltzmann's constant J/K
    k_eV = 8.617330350E-5  # Boltzmann's constant eV/K
    e = 1.602176620898E-19  # Coulomb

    delta = material_direct_allowed_gap - material_indirect_gap
    fr = exp(-delta / (k_eV * T))

    ## need to convert solar irradiance from Power/m**2(nm) into photon#/s*m**2(nm)
    ## power is Watt, which is Joule / s
    ## E = hc/wavelength
    ## at each wavelength, Power * (wavelength(m)/(h(Js)*c(m/s))) = ph#/s
    solar_spectra_photon_flux = solar_spectra_irradiance * (solar_spectra_wavelength_meters / (h * c))

    from scipy.integrate import simps
    ### Calculation of total solar power incoming
    P_in = simps(solar_spectra_irradiance, solar_spectra_wavelength)

    ### calculation of blackbody irradiance spectra
    ## units of W/(m**3), different than solar_spectra_irradiance!!! (This is intentional, it is for convenience)
    blackbody_irradiance = (2.0 * pi * h * c ** 2 / (solar_spectra_wavelength_meters ** 5)) * (
            1.0 / ((exp(h * c / (solar_spectra_wavelength_meters * k * T))) - 1.0))

    ## now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
    blackbody_photon_flux = blackbody_irradiance * (solar_spectra_wavelength_meters / (h * c))

    ## units of nm
    material_wavelength_for_absorbance_data = ((c * h_eV) / (material_eV_for_absorbance_data + 0.00000001)) * 10 ** 9

    ### absorbance interpolation onto each solar spectrum wavelength
    from scipy.interpolate import interp1d
    ## creates cubic spline interpolating function, set up to use end values as the guesses if leaving the region where data exists
    material_absorbance_data_function = interp1d(material_wavelength_for_absorbance_data, material_absorbance_data,
                                                 kind='cubic',
                                                 fill_value=(material_absorbance_data[0], material_absorbance_data[-1]),
                                                 bounds_error=False)

    material_interpolated_absorbance = zeros(len(solar_spectra_wavelength_meters))
    for i in range(0, len(solar_spectra_wavelength_meters)):
        ## Cutting off absorption data below the gap. This is done to deal with VASPs broadening of the calculated absorption data
        if solar_spectra_wavelength[i] < ((
                                                  c * h_eV) / material_direct_allowed_gap) * 10 ** 9 or cut_off_absorbance_below_direct_allowed_gap == False:
            material_interpolated_absorbance[i] = material_absorbance_data_function(solar_spectra_wavelength[i])

    absorbed_by_wavelength = 1.0 - exp(-2.0 * material_interpolated_absorbance * thickness)

    ###  Numerically integrating irradiance over wavelength array
    ## Note: elementary charge, not math e!  ## units of A/m**2   W/(V*m**2)
    J_0_r = e * pi * simps(blackbody_photon_flux * absorbed_by_wavelength, solar_spectra_wavelength_meters)

    J_0 = J_0_r / fr

    ###  Numerically integrating irradiance over wavelength array
    ## elementary charge, not math e!  ### units of A/m**2   W/(V*m**2)
    J_sc = e * simps(solar_spectra_photon_flux * absorbed_by_wavelength, solar_spectra_wavelength)

    #    J[i] = J_sc - J_0*(1 - exp( e*V[i]/(k*T) ) )   #This formula from the original paper has a typo!!
    #    J[i] = J_sc - J_0*(exp( e*V[i]/(k*T) ) - 1)    #Bercx chapter and papers have the correct formula (well, the correction on one paper)
    def J(V):
        J = J_sc - J_0 * (exp(e * V / (k * T)) - 1.0)
        return J

    def power(V):
        p = J(V) * V
        return p

    from scipy.optimize import minimize

    def neg_power(V):  # Minimizing the negative of power to optimize power
        return -power(V)

    results = minimize(neg_power, x0=np.array(
        [0.0]))  # passing a function as a variable, preconditioning at V=0, since this is a smooth function
    # results.x are the inputs as an array that give the highest value of Power
    V_Pmax = float(results.x)
    P_m = power(V_Pmax)

    efficiency = P_m / P_in

    ## This segment isn't needed for functionality at all, but can display a plot showing how the maximization of power by choosing the optimal voltage value works
    # if False:
    #     V = linspace(0, 2, 200)
    #     plt.figure
    #     plt.plot(V, J(V))
    #     plt.plot(V, power(V), linestyle='--')
    #     plt.show()
    #     print(power(V_Pmax))

    return efficiency


###########################################

if __name__ == '__main__':  # This condition is only true if you are executing this file. If you call this file as a library, this condition will be False

    #	material_direct_allowed_gap = 2.304
    #	material_indirect_gap = 2.294

    # This example file is sort of like CdTe's absorption data. (It's a line between the endpoints.)  Good enough for an example, but don't use it for actual scientific data.
    file_name = "absorption.dat"

    from numpy import loadtxt, arange, logspace

    ## Assumed data format has the energy (in eV) in the first column and the absorbance (in m^-1, not cm^-1!) in the 2nd column
    material_eV_for_absorbance_data, material_absorbance_data = loadtxt(file_name, usecols=[0, 1],
                                                                        unpack=True)  # absorbance in 2nd column

    SLME = calculate_SLME(material_eV_for_absorbance_data=material_eV_for_absorbance_data,
                          material_absorbance_data=material_absorbance_data,
                          material_direct_allowed_gap=material_direct_allowed_gap,
                          material_indirect_gap=material_indirect_gap,
                          thickness=5E-6, T=293.15)
    print('File :', file_name)
    print('Standard SLME :', SLME)

    # thickness_array = logspace(-8, -2, num=1000)
    # SLME_list = []
    # for thickness in thickness_array:
    #     SLME_list.append(calculate_SLME(material_eV_for_absorbance_data=material_eV_for_absorbance_data,
    #                                     material_absorbance_data=material_absorbance_data,
    #                                     material_direct_allowed_gap=material_direct_allowed_gap,
    #                                     material_indirect_gap=material_indirect_gap,
    #                                     thickness=thickness, T=293.15))

#	if True:
#		from matplotlib.pylab import *
#		figure(figsize = (3.2, 2.5), dpi = 200)
#		semilogx(thickness_array*1e6, SLME_list)
#		xlabel('Thickness ($\mu m$)')
#		ylabel('SLME')
#		minorticks_on()
#		ax=gca()
#		ax.set_xlim(0.01, ax.get_xlim()[1])
#		ax.set_ylim(0, ax.get_ylim()[1])
#		tight_layout(pad = 0.5)
#		show()
    #
    # np.savetxt('Thickness_list.txt', thickness_array * 1e6)
    # np.savetxt('SLME_list.txt', SLME_list)
    #
    # sat_thick = [10000.0] * 1000
    # k = 0
    # for i in range(0, 999):
    #     if SLME_list[i + 1] - SLME_list[i] < 0.0000005:
    #         sat_thick[k] = thickness_array[i] * 1e6
    #         k = k + 1
    # sat_thick_value = min(sat_thick[:])
    #
    # text_file = open("slme_data.txt", "w")
    # text_file.write("SLME at 5um thickness: %s" % SLME + '\n')
    # text_file.write("SLME at 1000um thickness: %s" % SLME_list[999] + '\n')
    # text_file.write("Saturated thickness value: %s" % sat_thick_value + '\n')
    # text_file.close()
    #
    # fig = plt.subplots(figsize=(8, 6), dpi=450)
    # plt.subplots_adjust(left=0.18, bottom=0.18, right=0.94, top=0.93)
    # plt.rc('font', family='Arial narrow')
    #
    # SLME_perc = [0.0] * 1000
    # for i in range(0, 1000):
    #     SLME_perc[i] = float(SLME_list[i]) * 100
    # plt.plot(thickness_array[:] * 1e6, SLME_perc[:], c='dodgerblue', lw=3, label='_nolegend_')
    #
    # plt.rc('xtick', labelsize=32)
    # plt.rc('ytick', labelsize=32)
    #
    # plt.xlabel('Thickness ($\mu$m)', c='k', fontname='Arial Narrow', size=32, labelpad=8)
    # plt.ylabel('SLME (%)', c='k', fontname='Arial Narrow', size=32, labelpad=12)
    #
    # plt.xscale('log')
    # # plt.xlim([0, 1000])
    # # plt.ylim([1.2, 2.2])
    # # plt.xticks([0, 25, 50, 75, 100])
    # # plt.xticks([0, 8, 16, 24, 32])
    # # plt.yticks([1.75, 2.00, 2.25, 2.50, 2.75])
    #
    # # plt.legend(loc='upper center', ncol=3, frameon=True, prop={'family':'Arial narrow','size':16})
    #
    # plt.savefig('plot_slme.pdf', dpi=450)
