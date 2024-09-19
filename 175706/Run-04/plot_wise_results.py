import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

from matplotlib.lines import Line2D

model = '1996FG3_binary_lowres'
model_diam = 1700.55
spectral_emissivity = 0.95

def wrap_data_steps(x):
    if x < 0:
        return x + 2925
    if x > 2925:
        return x - 2925
    return x


accepted_clones = pd.read_csv(f'{model}_accepted_clones.txt', delimiter='\t', names=['chi2', 'diameter', 'pv', 'therm_inertia', 'roughness', 'FCF'])
fitting_results = pd.read_csv(f'{model}_fitting_results.txt')
wavelengths_WISE = ['4.6028','11.0984','22.6405']
colours = ['r', 'g', 'b', 'k']

for run in [0, 1, 2, 3, 4]:
    if run == 0:
        wavs = [0, 1, 2]
    else:
        wavs = [0]

    for wavelength_id in wavs:

        obs_data = pd.read_csv(f'175706-Run-0{run}.txt', delimiter='\t', names=['LC_ID', 'Wavelength_ID', 'step', 'step_unc', 'flux', 'flux_unc'])
        obs_data.step += 1863
        obs_data['step'] = obs_data['step'].apply(wrap_data_steps)
        TI = round(float(fitting_results['Average Thermal Inertia'].values[0]))
        TI = math.ceil(TI / 10) * 10
        rough_flux_1 = pd.read_csv(f'{model}_rough_flux({TI}_{run}).txt', delimiter='\t', names=wavelengths_WISE, index_col=False)

        smooth_flux_1 = pd.read_csv(f'{model}_smooth_flux({TI}_{run}).txt', delimiter='\t', names=wavelengths_WISE, index_col=False)

        f_r = float(fitting_results['Average Roughness'].values[0])
        combined_flux_1 = f_r*rough_flux_1 + (1-f_r)*smooth_flux_1

        D = float(fitting_results['Average Diameter'].values[0])
        FCF = float(fitting_results['Average FCF'].values[0])
        scaled_flux_1 = combined_flux_1*np.power(D/model_diam, 2)*FCF
        scaled_flux_1['4.6028'] *= spectral_emissivity/0.95

        first = obs_data[(obs_data.LC_ID==0) & (obs_data.flux_unc != 1)]

        plt.rcParams.update({'font.size': 20})
        fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))
        ax1.errorbar(x=first[first.Wavelength_ID == wavelength_id].step/2925, y=first[first.Wavelength_ID == wavelength_id].flux/1e-16, yerr=first[first.Wavelength_ID == wavelength_id].flux_unc/1e-16, fmt='o', elinewidth=1, c=colours[wavelength_id], label=wavelengths_WISE[wavelength_id])
        ax1.plot(np.linspace(0, 1, 2925), scaled_flux_1[wavelengths_WISE[wavelength_id]][:2925]/1e-16, c=colours[wavelength_id])


        ax1.set_xlabel('Rotation Phase')
        ax1.set_ylabel(r'Flux [$10^{-16}\,Wm^{-2}\mu m^{-1}$]')
        ax1.set_xlim(0, 1)


        legend_elements = []
        for i, wavelength in enumerate(wavelengths_WISE):
            legend_elements.append(Line2D([0], [0], color=colours[i], lw=2, label=fr'${wavelength}\,\mu m$'))
            
        fig.legend(handles=legend_elements, loc='outside upper center', ncol=len(wavelengths_WISE), bbox_to_anchor = [0.5, 0])

        plt.savefig(f'{model}_Run-{run:02d}_Wav-{wavelength_id}_light_curves.png', bbox_inches='tight')
        plt.close()


# plt.rcParams.update({'font.size': 20})
# fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))
# for wavelength_id in [0,1,2]:
#     ax1.errorbar(x=first[first.Wavelength_ID == wavelength_id].step/2925, y=first[first.Wavelength_ID == wavelength_id].flux/1e-16, yerr=first[first.Wavelength_ID == wavelength_id].flux_unc/1e-16, fmt='o', elinewidth=1, c=colours[wavelength_id], label=wavelengths_WISE[wavelength_id])
#     ax1.plot(np.linspace(0, 1, 2925), scaled_flux_1[wavelengths_WISE[wavelength_id]][:2925]/1e-16, c=colours[wavelength_id])


# ax1.set_xlabel('Rotation Phase')
# ax1.set_ylabel(r'Flux [$10^{-16}\,Wm^{-2}\mu m^{-1}$]')
# ax1.set_xlim(0, 1)
# ax1.set_title('2010-04-30 (W2, W3, and W4 Fit)')
# ax1.set_yscale('log')


# legend_elements = []
# for i, wavelength in enumerate(wavelengths_WISE):
#     legend_elements.append(Line2D([0], [0], color=colours[i], lw=2, label=fr'${wavelength}\,\mu m$'))
    
# fig.legend(handles=legend_elements, loc='outside upper center', ncol=len(wavelengths_WISE), bbox_to_anchor = [0.5, 0])

# plt.savefig(f'{model}_log_light_curves.png', bbox_inches='tight')
# plt.close()

all_clones = pd.read_csv(f'{model}_all_clones.txt', delimiter='\t', names=['chi2', 'diameter', 'pv', 'therm_inertia', 'roughness', 'FCF'])
accepted_clones = pd.read_csv(f'{model}_accepted_clones.txt', delimiter='\t', names=['chi2', 'diameter', 'pv', 'therm_inertia', 'roughness', 'FCF'])
one_sigma_clones = all_clones[all_clones.chi2 < (all_clones.chi2.min() + 3.53)]

diameters = sorted(all_clones.diameter.unique())
thermal_inertias = sorted(all_clones.therm_inertia.unique())
roughnesses = sorted(all_clones.roughness.unique())

chi2_plane = np.zeros(shape=(len(thermal_inertias), len(roughnesses)))

for i, thermal_inertia in enumerate(thermal_inertias):
    for j, roughness in enumerate(roughnesses):
        clones = all_clones[(all_clones.therm_inertia == thermal_inertia) & (all_clones.roughness == roughness)]
        chi2_plane[i, j] = clones.chi2.min()

xx, yy = np.meshgrid(thermal_inertias, roughnesses)
levels = [all_clones.chi2.min() + 3.53, all_clones.chi2.min() + 8.02, all_clones.chi2.min() + 14.16]
plt.rcParams.update({'font.size': 20})
plt.figure(1, (10, 8))
CS = plt.contour(xx, yy, chi2_plane.T, levels=levels)
fmt = {}
strs = [r'$1\sigma$', r'$2\sigma$', r'$3\sigma$']
for l, s in zip(CS.levels, strs):
    fmt[l] = s
plt.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=10)
plt.xlabel(r'Thermal Inertia [J$\,$m$^{-2}\,$K$^{-1}\,$s$^{-1/2}$]')
plt.ylabel('Roughness Fraction')
plt.savefig('Thermal_inertia-Roughness-contours.png', bbox_inches='tight')
plt.close()

print(fr"$\text{{Diameter}} = {accepted_clones.diameter.mean():.0f} \pm {accepted_clones.diameter.std():.0f}\,\text{{km}}$")
print(fr"$\text{{Thermal Inertia}} = {accepted_clones.therm_inertia.mean():.0f} \pm {accepted_clones.therm_inertia.std():.0f}\,\text{{J}}\,\text{{m}}^{{-2}}\,\text{{K}}^{{-1}}\,\text{{s}}^{{-1/2}}$")
print(fr"$\text{{Roughness Fraction}} = {accepted_clones.roughness.mean():.2f} \pm {accepted_clones.roughness.std():.2f}$")
