
import os, shutil, sys
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def plot_pinwise_errors_BWR_assembly(errors_rates, assembly_model, assembly_id, name_compo, calculation_opt, fig_name, evaluation, cmap="Spectral"):
    """
    plot_pinwise_errors_BWR_assembly: Plot pin-wise relative differences in rates between Dragon5 and Serpent2 for a BWR assembly.
    The function takes in a 2D array of relative differences (errors_rates) for each pin and energy group, the assembly geometry (assembly_map), and saves a heatmap of the relative differences for each energy group.
    Parameters:
    - errors_rates: A 2D numpy array of shape (ngroups, ncells)
    - assembly_model: A CartesianAssemblyModel object containing the geometry and material information of the assembly.
    - assembly_id: A string identifier for the assembly (e.g., "GE14_DOM")
    - name_compo: name of the DRAGON compo file parsed
    - calculation_opt: name of the calculation option (e.g., "RSE_IC_1L_MOC")
    - fig_name: name of the quantity plotted in the figure (e.g., "Fiss_over_abs_diff" or "Fiss_rates")
    - evaluation: name of the nuclear data evaluation used (e.g., "endfb8r1")
    
    """
    a = os.path.exists(f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{name_compo}")
    print(f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{name_compo} already exists: {a}")
    if not a:
        os.makedirs(f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{name_compo}")
        
    results_dir = f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{name_compo}"
        
    # Get lattice info from assembly model
    
    lattice_info = assembly_model.get_postprocessing_lattice_info()
    ordered_pin_indices = lattice_info['ordered_pin_indices']
    pin_idx_on_axis = lattice_info['pin_idx_on_symmetry_axis']
    n_unique_pins = lattice_info['n_unique_pins']
    lattice = assembly_model.lattice
    lattice_description = assembly_model.lattice_description
    
    print(f"Assembly {assembly_id} has {n_unique_pins} unique pins, with indices {ordered_pin_indices}, and pin indices on symmetry axis {pin_idx_on_axis}.")
    print(f"Lattice shape: {np.array(lattice).shape}, with values:\n{lattice}")
    print(f"Lattice description: {lattice_description}")
    ngroups = errors_rates.shape[0]
    ncells = errors_rates.shape[1]
    for gr in range(ngroups):
        
        error_grid = np.full_like(lattice_description, np.nan, dtype=float)
        k=0
        for i in range(error_grid.shape[0]):
            for j in range(error_grid.shape[1]):
                print(error_grid.shape)
                print(f"Processing pin at position ({i}, {j}) in the lattice...")
                ROD_ID = lattice_description[i][j]
                if ROD_ID in assembly_model.non_fuel_rod_ids:
                    continue
                    
                if assembly_model.check_diagonal_symmetry() == "anti-diagonal":
                    # Below anti-diagonal
                    if i + j <= error_grid.shape[0] - 1:
                        error_grid[i, j] = errors_rates[gr,k]
                        # Mirror if not on symmetry line
                        if i + j < error_grid.shape[0] - 1:
                            ii = error_grid.shape[0] - 1 - j
                            jj = error_grid.shape[0] - 1 - i
                            error_grid[ii, jj] = errors_rates[gr,k]
                        k += 1
                elif assembly_model.check_diagonal_symmetry() == "main-diagonal":
                    # Below main diagonal
                    if i >= j:
                        error_grid[i, j] = errors_rates[gr,k]
                        # Mirror if not on symmetry line
                        if i > j:
                            ii = j
                            jj = i
                            error_grid[ii, jj] = errors_rates[gr,k]
                        k += 1
        # Flip once
        plot_grid = error_grid #np.flipud(error_grid)

        plt.figure(figsize=(10,10))

        im = plt.imshow(
            plot_grid,
            origin="lower",
            cmap=cmap,
            vmin=np.min(errors_rates[gr]),
            vmax=np.max(errors_rates[gr])
        )

        plt.colorbar(im, label="Relative Difference (%)")

        nrows, ncols = plot_grid.shape
            
        for i in range(nrows):
            for j in range(ncols):
                val = plot_grid[i, j]
                if fig_name == "fiss_over_abs_diff":
                    formatted_val = f"{val:.2f}"
                elif fig_name == "fission_rates_diff":
                    formatted_val = f"{val:.2f}%"
                if not np.isnan(val):
                    plt.text(
                        j, i,
                        formatted_val,
                        ha="center",
                        va="center",
                        color="black",
                        fontsize=12,
                        fontweight="bold"
                    )

        # Grid lines
        ax = plt.gca()
        ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
        ax.grid(which="minor", color="black", linestyle="-", linewidth=0.8)

        ax.tick_params(which="both", bottom=False, left=False,
               labelbottom=False, labelleft=False)

        # Clip the grid strictly to the data area
        ax.set_xlim(-0.5, ncols-0.5)
        ax.set_ylim(-0.5, nrows-0.5)
                    
        if gr == 0:
            group_id = "thermal"
        elif gr == 1:
            group_id = "fast"
        if fig_name == "fiss_over_abs_diff":
            plt.title(f'(D5-S2) Relative differences on $\\tau_f/\\tau_a$, {group_id} group', fontsize=16)
        elif fig_name == "fission_rates_diff":
            plt.title(f'(D5-S2) Relative differences on fission rates, {group_id} group', fontsize=16)
        # Show the plot
        plt.tight_layout()
        plt.savefig(f"{results_dir}/assembly_map_{fig_name}_g{gr+1}.png", dpi=300)
        plt.close()
        

    return


def plot_spectrum_comparison(energy_mesh, FLUX_groups, S2_spectra, assembly_id, CPO_name, calculation_opt, evaluation):
    n_groups = len(FLUX_groups)
    try:
        S2_spectrum = S2_spectra[f"{n_groups}g"]
    except KeyError:
        print(f"Error: Serpent2 spectrum for {n_groups} groups not found in S2_spectra dictionary.")
        return    
    
    a = os.path.exists(f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{CPO_name}")
    print(f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{CPO_name} already exists: {a}")
    if not a:
        os.makedirs(f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{CPO_name}")
    results_dir = f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{CPO_name}"
    
    # Compute lethargy mesh for plotting spectrum on lethargy scale
    lethargy = np.log(np.max(energy_mesh) / energy_mesh)
    # Normalize DRAGON spectrum to Serpent2 spectrum for better comparison
    normalization_factor = np.sum(S2_spectrum) / np.sum(FLUX_groups)
    FLUX_groups = FLUX_groups * normalization_factor
    # plot flux spectrum
    plt.figure(figsize=(10, 6))
    # FLUX has 295 values = 1 per group, extend so that it is ploted as piecewise constant between energy groups
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(FLUX_groups, 2), where='post', label='DRAGON5 Flux Spectrum', color='blue', linewidth=2)
    plt.plot(np.repeat(energy_mesh, 2)[1:-1], np.repeat(S2_spectrum, 2), label='Serpent2 Flux Spectrum', color='orange', linewidth=2)
    # plot vertical line at 0.625 eV
    #plt.vlines([0.625, 10, 1000], 0.0, 0.032, colors = ["red"])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Flux (n/cm²/s)')
    plt.title('Flux Spectrum Comparison')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{results_dir}/flux_spectrum_comparison.png")
    plt.close()
    
    # plot flux spectrum on lethargy scale
    plt.figure(figsize=(10, 6))
    # FLUX has 295 values = 1 per group, extend so that it is ploted as piecewise constant between energy groups
    plt.step(np.repeat(lethargy, 2)[1:-1], np.repeat(FLUX_groups, 2), where='post', label='DRAGON5 Flux Spectrum', color='blue', linewidth=2)
    plt.plot(np.repeat(lethargy, 2)[1:-1], np.repeat(S2_spectrum, 2), label='Serpent2 Flux Spectrum', color='orange', linewidth=2)
    # plot vertical line at 0.625 eV
    #plt.vlines([0.625, 10, 1000], 0.0, 0.032, colors = ["red"])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Lethargy')
    plt.ylabel('Flux (n/cm²/s)')
    plt.title('Flux Spectrum Comparison')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{results_dir}/flux_spectrum_comparison_lethargy.png")
    plt.close()

    delta_flx = FLUX_groups - S2_spectrum
    absolute_diff = np.abs(delta_flx)
    delta_rel_flx = (FLUX_groups - S2_spectrum) * 100 / S2_spectrum

    # Plot relative difference in flux
    plt.figure(figsize=(10, 6))
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(delta_rel_flx, 2), where='post', label='Relative Difference in Flux', color='green', linewidth=2)
    plt.xscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Relative Difference (%)')
    plt.title('Relative Difference in Flux Spectrum')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{results_dir}/relative_difference_flux_spectrum.png")
    plt.close()

    # Plot absolute difference in flux
    plt.figure(figsize=(10, 6))
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(absolute_diff, 2), where='post', label='Absolute Difference in Flux', color='red', linewidth=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Flux Absolute Difference')
    plt.title('Absolute Difference in Flux Spectrum')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{results_dir}/absolute_difference_flux_spectrum.png")
    plt.close()

    # Plot difference in flux
    plt.figure(figsize=(10, 6))
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(delta_flx, 2), where='post', label='Difference in Flux', color='red', linewidth=2)
    plt.xscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Flux Difference')
    plt.title('Difference in Flux Spectrum')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{results_dir}/difference_flux_spectrum.png")
    plt.close()

    
def plot_U238_abs_rates(energy_mesh, D5_U238_abs_rates, S2_U238_abs_rates_dict, assembly_id, CPO_name, calculation_opt, evaluation):
    """ 
    plot spatially integrated U238 absorption rates on energy mesh
    """
    results_dir = f"DRAGON_starterDD_RATES_{assembly_id}_{evaluation}/{calculation_opt}/{CPO_name}"
    a = os.path.exists(results_dir)
    print(f"{results_dir} already exists: {a}")
    if not a:
        os.makedirs(results_dir)
    
    n_groups = len(D5_U238_abs_rates)
    # recover U238 absorption rates collapsed to n groups from S2_U238_abs_rates dictionary
    cropping_window = [295-206, 295-52] 
    try:        
        S2_U238_abs_rates = S2_U238_abs_rates_dict[f"{n_groups}g"]
    except KeyError:
        print(f"Error: Serpent2 U238 absorption rates for {n_groups} groups not found in S2_U238_abs_rates_dict dictionary.")
        return
    # Normalize DRAGON rates to Serpent2 rates for better comparison
    normalization_factor = np.sum(S2_U238_abs_rates) / np.sum(D5_U238_abs_rates)
    D5_U238_abs_rates = D5_U238_abs_rates * normalization_factor
    diff_abs_rates = D5_U238_abs_rates - S2_U238_abs_rates
    rel_diff_abs_rates = (D5_U238_abs_rates - S2_U238_abs_rates) * 100 / S2_U238_abs_rates
    max_rel_diff = np.max(np.abs(rel_diff_abs_rates))
    avg_rel_diff = np.mean(np.abs(rel_diff_abs_rates))
    rms_rel_diff = np.sqrt(np.mean(rel_diff_abs_rates**2))
    # group where max relative difference is observed
    group_max_rel_diff = n_groups - np.argmax(np.abs(rel_diff_abs_rates))
    print(f"Max relative difference in U238 absorption rates: {max_rel_diff:.2f}%, in group {group_max_rel_diff}")
    print(f"Average relative difference in U238 absorption rates: {avg_rel_diff:.2f}%")
    print(f"RMS relative difference in U238 absorption rates: {rms_rel_diff:.2f}%")
    # plot absorption rates
    plt.figure(figsize=(10, 6))
    # D5_U238_abs_rates has 295 values = 1 per group, extend so that it is ploted as piecewise constant between energy groups
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(D5_U238_abs_rates, 2), where='post', label='DRAGON5 U238 Absorption Rates', color='blue', linewidth=2)
    plt.plot(np.repeat(energy_mesh, 2)[1:-1], np.repeat(S2_U238_abs_rates, 2), label='Serpent2 U238 Absorption Rates', color='orange', linewidth=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Absorption Rate (1/s)')
    plt.title('U238 Absorption Rates Comparison')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{results_dir}/U238_absorption_rates_comparison.png", dpi=300)
    plt.close()
    
    # plot relative difference in absorption rates
    plt.figure(figsize=(10, 6))
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(rel_diff_abs_rates, 2), where='post', label='Relative Difference in U238 Absorption Rates', color='green', linewidth=2)
    plt.xscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Relative Difference (%)')
    plt.title('Relative Difference in U238 Absorption Rates')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{results_dir}/relative_difference_U238_absorption_rates.png", dpi=300)
    plt.close()
    
    if n_groups == 295:
        cropping_window = [295-206, 295-52] 
        x = energy_mesh[cropping_window[0]:cropping_window[1]]
        y = rel_diff_abs_rates[cropping_window[0]:cropping_window[1]]
        max_rel_diff_zoom = np.max(np.abs(y))
        avg_rel_diff_zoom = np.mean(np.abs(y))
        rms_rel_diff_zoom = np.sqrt(np.mean(y**2))
        group_max_rel_diff_zoom = n_groups - cropping_window[0] - np.argmax(np.abs(y))
        print(f"Zoom on energy range of interest [{energy_mesh[cropping_window[1]]:.2e} eV, {energy_mesh[cropping_window[0]]:.2e} eV]:")
        print(f"Max relative difference in U238 absorption rates: {max_rel_diff_zoom:.2f}%, in group {group_max_rel_diff_zoom}, from E1: {energy_mesh[group_max_rel_diff_zoom-1]:.2e} eV to E2: {energy_mesh[group_max_rel_diff_zoom]:.2e} eV")
        print(f"Average relative difference in U238 absorption rates: {avg_rel_diff_zoom:.2f}%")
        print(f"RMS relative difference in U238 absorption rates: {rms_rel_diff_zoom:.2f}%")

        plt.figure(figsize=(10, 6))
        plt.step(x, y, where='post', linewidth=2)

        plt.xscale('log')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Relative Difference (%)')
        plt.title('Relative Difference in U238 Absorption Rates (g52 to g206)')
        plt.grid()
        plt.tight_layout()
        plt.savefig(f"{results_dir}/relative_difference_U238_absorption_rates_zoom.png", dpi=300)
        plt.close()