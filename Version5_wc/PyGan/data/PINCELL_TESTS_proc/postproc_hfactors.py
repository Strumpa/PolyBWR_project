## Post treat H-FACTOR vs ENERGY for U5, U8 and Pu9.
# Investigate H-FACTOR records in :
# 1) Microlib prior to self-shielding,
# 2) self-shielded 
import matplotlib.pyplot as plt
import numpy as np

def get_energy_midpoints(energy_bounds):
    """Calculate midpoint energies for each group."""
    midpoints = np.sqrt(energy_bounds[:-1] * energy_bounds[1:])
    return midpoints

def plot_microlib_info(pyLib, save_name, save_dir):
    """
    pyLib : lcm object containing mirolib information from LIB: module 
    """

    energy_mesh = pyLib["ENERGY"]
    h_factors = {}
    nftot = {}
    nTOT = {}
    nWT0 = {}

    isos = ["U235", "U238"]#["U234", "U235", "U238", "Pu239"]
    
    for iso_idx in range(len(pyLib["ISOTOPESDENS"])):
        density = pyLib["ISOTOPESDENS"][iso_idx]
        if density > 0.0:
            print(f"For isotope with index : {iso_idx}, density = {pyLib["ISOTOPESDENS"][iso_idx]}")
            print(f"H-FACTOR = {pyLib["ISOTOPESLIST"][iso_idx]["H-FACTOR"]}")
            print(f"ISOTOPESLIST[iso_idx]['ALIAS'][0:5].strip() = {pyLib["ISOTOPESLIST"][iso_idx]['ALIAS'][0:5].strip()}")
            alias = pyLib["ISOTOPESLIST"][iso_idx]['ALIAS'][0:5].strip()
            if alias in isos:
                h_factors[alias] = pyLib["ISOTOPESLIST"][iso_idx]["H-FACTOR"]
                nftot[alias] = pyLib["ISOTOPESLIST"][iso_idx]["NFTOT"]
                nTOT[alias] = pyLib["ISOTOPESLIST"][iso_idx]["NTOT0"]
                nWT0[alias] = pyLib["ISOTOPESLIST"][iso_idx]["NWT0"]

    for iso_name in isos:
        if iso_name in h_factors.keys():

            fig, axes = plt.subplots(3, 1, figsize=(12, 10))
            fig.suptitle(f'{iso_name} H-FACTOR / NFTOT data Recovered from LIB - {save_name}',
                        fontsize=14, fontweight='bold')

            # Plot H-FACTOR
            energy_midpoints = get_energy_midpoints(energy_mesh)
            h_factors_iso = h_factors[iso_name]
            axes[0].loglog(energy_midpoints, h_factors_iso, 'b-o', markersize=4,
                        linewidth=1.5, label='H-FACTOR')
            axes[0].set_xlabel('Energy (eV)', fontsize=11)
            axes[0].set_ylabel('H-FACTOR (eV·b)', fontsize=11)
            axes[0].grid(True, alpha=0.3)
            axes[0].legend()
               
            # Plot NFTOT
            nftot_iso = nftot[iso_name]
            axes[1].loglog(energy_midpoints, nftot_iso, 'r-s', markersize=4,
                        linewidth=1.5, label='NFTOT')
            axes[1].set_xlabel('Energy (eV)', fontsize=11)
            axes[1].set_ylabel('NFTOT (b)', fontsize=11)
            axes[1].grid(True, alpha=0.3, which='both')
            axes[1].legend()

            # Plot RATIO
            ratios = h_factors_iso / nftot_iso
            #axes[2].semilogx(energy_midpoints, ratios, 'g-^', markersize=4,
            #                linewidth=1.5, label='H-FACT/NFTOT Ratio')
            axes[2].loglog(energy_midpoints, ratios, 'g-^', markersize=4,
                            linewidth=1.5, label='H-FACT/NFTOT Ratio')
            axes[2].set_xlabel('Energy (eV)', fontsize=11)
            axes[2].set_ylabel('Ratio (eV)', fontsize=11)
            axes[2].grid(True, alpha=0.3)
            axes[2].legend()

            plt.tight_layout()
            output_file = f'{save_dir}/hfact_{save_name}_{iso_name}.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Plot saved as '{output_file}'")
        



def plot_self_shielded_microlib_info(pyLIB2, save_name, save_dir):
    """
    pyLIB2 : lcm object containing mirolib information from LIB: module 
    """

    
    print(pyLIB2["STATE-VECTOR"])
    energy_mesh = pyLIB2["ENERGY"]
    print(energy_mesh)
    print(len(pyLIB2["ISOTOPESLIST"]))
    h_factors = {}
    nftot = {}
    nTOT={}
    isos = ["U234", "U235", "U238", "Pu239"]
    
    #print(pyLIB2["ISOTOPESUSED"])
    print(f"isotopes list 0 : {pyLIB2["ISOTOPESLIST"][0]}")
    #print(pyLIB2["ISOTOPESDENS"][0])
    for iso_idx in range(len(pyLIB2["ISOTOPESDENS"])):
        #print(f"Isotope index = {iso_idx}, alias = {pyLib["ISOTOPESLIST"][iso_idx]["ALIAS"]}")
        #isotope = pyLib["ISOTOPESLIST"][iso_idx]["ALIAS"][0:5].strip()
        density = pyLIB2["ISOTOPESDENS"][iso_idx]
        if density > 0.0:
            print(f"For isotope with index : {iso_idx}, density = {pyLIB2["ISOTOPESDENS"][iso_idx]}")
            print(f"H-FACTOR = {pyLIB2["ISOTOPESLIST"][iso_idx]["H-FACTOR"]}")
            print(f"ISOTOPESLIST[iso_idx]['ALIAS'][0:5].strip() = {pyLIB2["ISOTOPESLIST"][iso_idx]['ALIAS'][0:5].strip()}")
            alias = pyLIB2["ISOTOPESLIST"][iso_idx]['ALIAS'][0:5].strip()
            if alias in isos:
                h_factors[alias] = pyLIB2["ISOTOPESLIST"][iso_idx]["H-FACTOR"]
                nftot[alias] = pyLIB2["ISOTOPESLIST"][iso_idx]["NFTOT"]
                nTOT[alias] = pyLIB2["ISOTOPESLIST"][iso_idx]["NTOT0"]

    for iso_name in isos:
        if iso_name in h_factors.keys():

            fig, axes = plt.subplots(3, 1, figsize=(12, 10))
            fig.suptitle(f'{iso_name} Data Recovered from self-shielded LIB - {save_name}',
                        fontsize=14, fontweight='bold')

            # Plot H-FACTOR
            energy_midpoints = get_energy_midpoints(energy_mesh)
            h_factors_iso = h_factors[iso_name]
            axes[0].loglog(energy_midpoints, h_factors_iso, 'b-o', markersize=4,
                        linewidth=1.5, label='H-FACTOR')
            axes[0].set_xlabel('Energy (eV)', fontsize=11)
            axes[0].set_ylabel('H-FACTOR (eV·b)', fontsize=11)
            axes[0].grid(True, alpha=0.3)
            axes[0].legend()
               
            # Plot NFTOT
            nftot_iso = nftot[iso_name]
            axes[1].loglog(energy_midpoints, nftot_iso, 'r-s', markersize=4,
                        linewidth=1.5, label='NFTOT')
            axes[1].set_xlabel('Energy (eV)', fontsize=11)
            axes[1].set_ylabel('NFTOT (b)', fontsize=11)
            axes[1].grid(True, alpha=0.3, which='both')
            axes[1].legend()

            # Plot RATIO
            ratios = h_factors_iso / nftot_iso
            print(ratios)
            #axes[2].semilogx(energy_midpoints, ratios, 'g-^', markersize=4,
            #                linewidth=1.5, label='H-FACT/NFTOT Ratio')
            axes[2].loglog(energy_midpoints, ratios, 'g-^', markersize=4,
                            linewidth=1.5, label='H-FACT/NFTOT Ratio')
            axes[2].set_xlabel('Energy (eV)', fontsize=11)
            axes[2].set_ylabel('Ratio (eV)', fontsize=11)
            axes[2].grid(True, alpha=0.3)
            axes[2].legend()

            plt.tight_layout()
            output_file = f'{save_dir}/self_shielded_hfact_{save_name}_{iso_name}.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Plot saved as '{output_file}'")
