# Python3 script grouping plotters for OECD NEA Phase IIIB benchmark calculations.
import re
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_spectrum_comparison(energy_mesh, FLUX_295groups, S2_spectrum, cpo_name, results_dir):
    results_dir = f"{results_dir}/{cpo_name}"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # plot flux spectrum
    plt.figure(figsize=(10, 6))
    # FLUX has 295 values = 1 per group, extend so that it is ploted as piecewise constant between energy groups
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(FLUX_295groups, 2), where='post', label='DRAGON5 Flux Spectrum', color='blue', linewidth=2)
    plt.plot(np.repeat(energy_mesh, 2)[1:-1], np.repeat(S2_spectrum, 2), label='Serpent2 Flux Spectrum', color='orange', linewidth=2)
    # plot vertical line at 0.625 eV
    plt.vlines([0.625, 10, 1000], 0.0, 0.032, colors = ["red"])
    plt.xscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Flux (n/cm²/s)')
    plt.title('Flux Spectrum Comparison')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{results_dir}/flux_spectrum_comparison.png')
    plt.close()

    delta_flx = FLUX_295groups - S2_spectrum
    absolute_diff = np.abs(delta_flx)
    delta_rel_flx = (FLUX_295groups - S2_spectrum) * 100 / S2_spectrum

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
    plt.savefig(f'{results_dir}/relative_difference_flux_spectrum.png')
    plt.close()

    # Plot absolute difference in flux
    plt.figure(figsize=(10, 6))
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(absolute_diff, 2), where='post', label='Absolute Difference in Flux', color='red', linewidth=2)
    plt.xscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Flux Absolute Difference (n/cm²/s)')
    plt.title('Absolute Difference in Flux Spectrum')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{results_dir}/absolute_difference_flux_spectrum.png')
    plt.close()

    # Plot difference in flux
    plt.figure(figsize=(10, 6))
    plt.step(np.repeat(energy_mesh, 2)[1:-1], np.repeat(delta_flx, 2), where='post', label='Difference in Flux', color='red', linewidth=2)
    plt.xscale('log')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Flux Difference (n/cm²/s)')
    plt.title('Difference in Flux Spectrum')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{results_dir}/difference_flux_spectrum.png')
    plt.close()

def plot_eighth_lattice(compo_name, delta_reaction_rates, results_dir):
    """
    Plot the fission rates for the eighth lattice from Serpent2 and DRAGON5.
    """
    
    a = os.path.exists(f"{results_dir}/{compo_name}")
    if not a:
        os.makedirs(f"{results_dir}/{compo_name}")
    # Create a figure and axis
    fig1, ax1 = plt.subplots(figsize=(10, 10))
    fig2, ax2 = plt.subplots(figsize=(10, 10))
    cmap = "Spectral"
    color_map = plt.cm.get_cmap(cmap)
    # group 1 norm
    min_error1 = min(list(delta_reaction_rates[0]))
    max_error1 = max(list(delta_reaction_rates[0]))
    norm1 = plt.Normalize(min_error1, max_error1)
    # group 2 norm
    min_error2 = min(list(delta_reaction_rates[1]))
    max_error2 = max(list(delta_reaction_rates[1]))
    norm2 = plt.Normalize(min_error2, max_error2)

    # Define the grid size
    grid_size = 4


    # Define the colors for each cell
    colors = {
        "1": "blue",
        "2": "green",
        "3": "yellow",
        "4": "orange",
        "5": "red",
        "6": "purple",
        "7": "pink",
        "8": "cyan",
        "9" : "white",
        "WR": "black",
        "sy": "lightgrey"
    }
    assembly_map = np.array([
        ["sy", "sy", "sy", "1"],
        ["sy", "sy", "5", "2"],
        ["sy", "8", "6", "3"],
        ["WR", "9", "7", "4"],
        ])

    # Loop through the assembly map and plot each cell
    for row in range(assembly_map.shape[0]):
        for col in range(assembly_map.shape[1]):
            cell_name = assembly_map[row][col]
            if cell_name in colors.keys():
                color = colors[cell_name]
                # Get the reaction rate for the current cell
                
                # Get the reaction rate for the current cell
                # Plot the cell with the corresponding color
                if cell_name == "sy" or cell_name == "WR":
                    print("skip") 
                else:
                    cell_name_idx = int(cell_name) - 1  # Convert cell name to index (0-based)
                    delta_rate_g1 = delta_reaction_rates[0][cell_name_idx]
                    delta_rate_g2 = delta_reaction_rates[1][cell_name_idx]
                    color1 = color_map(norm1(delta_rate_g1))  # Get color from colormap
                    color2 = color_map(norm2(delta_rate_g2))
                    rect1 = plt.Rectangle((col, grid_size - row - 1), 1, 1, color=color1, edgecolor='black')
                    rect2 = plt.Rectangle((col, grid_size - row - 1), 1, 1, color=color2, edgecolor='black')
                    ax1.add_patch(rect1)
                    ax2.add_patch(rect2)
                    # Add the reaction rate text inside the cell
                    print(f"cell_name: {cell_name}, delta_rate_g1: {delta_rate_g1}, rate_g2: {delta_rate_g2}")
                    ax1.text(col + 0.5, grid_size - row - 0.5, f"cell pos = {cell_name} {delta_rate_g1:.2f}", ha='center', va='center', fontsize=12, color='black')
                    ax2.text(col + 0.5, grid_size - row - 0.5, f"cell pos = {cell_name} {delta_rate_g2:.2f}", ha='center', va='center', fontsize=12, color='black')
    # Set the limits and aspect ratio
    ax1.set_xlim(0, grid_size)
    ax1.set_ylim(0, grid_size)
    ax1.set_aspect('equal')
    # Set the title and labels
    ax1.set_title(f"Relative difference on U235 fission rates, thermal group", fontsize=16)
    ax1.set_xlabel("X-axis", fontsize=12)
    ax1.set_ylabel("Y-axis", fontsize=12)
    # Hide the axes
    ax1.set_xticks([])
    ax1.set_yticks([])
    # Show the plot
    plt.tight_layout()
    fig1.savefig(f"{results_dir}/{compo_name}/OECD_NEA_PHASE_IIIB_fiss_rates_assembly_map_g1.png", dpi=300)
    plt.close(fig1)

    ax2.set_xlim(0, grid_size)
    ax2.set_ylim(0, grid_size)
    ax2.set_aspect('equal')
    # Set the title and labels
    ax2.set_title(f"Relative difference on U235 fission rates, fast group", fontsize=16)
    ax2.set_xlabel("X-axis", fontsize=12)
    ax2.set_ylabel("Y-axis", fontsize=12)
    # Hide the axes
    ax2.set_xticks([])
    ax2.set_yticks([])
    # Show the plot
    plt.tight_layout()
    fig2.savefig(f"{results_dir}/{compo_name}/OECD_NEA_PHASE_IIIB_fiss_rates_assembly_map_g2.png", dpi=300)
    plt.close(fig2)
    return