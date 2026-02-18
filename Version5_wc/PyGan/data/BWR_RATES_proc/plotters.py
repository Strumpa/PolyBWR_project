
import os, shutil, sys
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def plot_BWR_assembly(reaction_rates, name_case, name_compo, fig_name, unfold_symmetry):
        """
        Assembly map for ATRIUM-10

        C1 C2 C3 C5 C6 C5 C4 C3 C2 C1
        C2 C4 C7 C6 C7 C6 C6 C7 C4 C2
        C3 C7 C6 C6 C6 C7 C6 C6 C7 C3
        C5 C6 C6 C6 C6 C6 C7 C6 C5 C4
        C6 C7 C6 C6 WH WH WH C4 C6 C4
        C5 C6 C7 C6 WH WH WH C3 C7 C4
        C4 C6 C6 C7 WH WH WH C4 C4 C4
        C3 C7 C6 C6 C4 C3 C4 C4 C8 C3
        C2 C4 C7 C5 C6 C7 C4 C8 C4 C2
        C1 C2 C3 C4 C4 C4 C4 C3 C2 C1

        WH = Water Hole
        Cn = Fuel Cell number n

        Parameters:
        - reaction_rates (dict): Dictionary containing reaction rates for each cell.
        - name_case (str): Name of the geometry.
        - name_compo (str): Name of the MULTICOMPO obect with calculation parameter identifiers.
        - fig_name (str): Name of the figure to be saved.
        - unfold_symmetry (bool): If True, apply symmetry to the assembly map.
        """
        a = os.path.exists(f"DRAGON_RATES_{name_case}/{name_compo}")
        if not a:
            os.makedirs(f"DRAGON_RATES_{name_case}/{name_compo}")
        # Create a figure and axis
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        fig2, ax2 = plt.subplots(figsize=(10, 10))

        # Define the grid size
        grid_size = 10


        # Define the colors for each cell
        colors = {
            "C1": "blue",
            "C2": "green",
            "C3": "yellow",
            "C4": "orange",
            "C5": "red",
            "C6": "purple",
            "C7": "pink",
            "C8": "cyan",
            "WH": "white",
            "sy": "lightgrey"
        }
        if unfold_symmetry:
            # Define the assembly map
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["C2", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["C3", "C7", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["C5", "C6", "C6", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["C6", "C7", "C6", "C6", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["C5", "C6", "C7", "C6", "WH", "WH", "WH", "C3", "C7", "C4"],
                ["C4", "C6", "C6", "C7", "WH", "WH", "WH", "C4", "C4", "C4"],
                ["C3", "C7", "C6", "C6", "C4", "C3", "C4", "C4", "C8", "C3"],
                ["C2", "C4", "C7", "C5", "C6", "C7", "C4", "C8", "C4", "C2"],
                ["C1", "C2", "C3", "C4", "C4", "C4", "C4", "C3", "C2", "C1"]
            ])
        else:
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["sy", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["sy", "sy", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["sy", "sy", "sy", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["sy", "sy", "sy", "sy", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "WH", "WH", "C3", "C7", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "WH", "C4", "C4", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C8", "C3"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C2"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C1"]
                ])

        # Loop through the assembly map and plot each cell
        for row in range(assembly_map.shape[0]):
            for col in range(assembly_map.shape[1]):
                cell_name = assembly_map[row][col]
                if cell_name in colors.keys():
                    color = colors[cell_name]
                    # Get the reaction rate for the current cell
                    if cell_name in reaction_rates.keys():
                        rate_g1 = reaction_rates[cell_name][0]
                        rate_g2 = reaction_rates[cell_name][1]
                    else:
                        rate_g1 = cell_name
                        rate_g2 = cell_name
                    # Get the reaction rate for the current cell
                    # Plot the cell with the corresponding color
                    if cell_name == "sy" or cell_name == "WH":
                        print("skip") 
                    else:
                        rect1 = plt.Rectangle((col, grid_size - row - 1), 1, 1, color=color, edgecolor='black')
                        rect2 = plt.Rectangle((col, grid_size - row - 1), 1, 1, color=color, edgecolor='black')
                        ax1.add_patch(rect1)
                        ax2.add_patch(rect2)
                        # Add the reaction rate text inside the cell
                        ax1.text(col + 0.5, grid_size - row - 0.5, f"{rate_g1:.2f}", ha='center', va='center', fontsize=8, color='black')
                        ax2.text(col + 0.5, grid_size - row - 0.5, f"{rate_g2:.2f}", ha='center', va='center', fontsize=8, color='black')
        # Set the limits and aspect ratio
        ax1.set_xlim(0, grid_size)
        ax1.set_ylim(0, grid_size)
        ax1.set_aspect('equal')
        # Set the title and labels
        ax1.set_title(f"Assembly Map for {name_case} - {name_compo}, fast group", fontsize=16)
        ax1.set_xlabel("X-axis", fontsize=12)
        ax1.set_ylabel("Y-axis", fontsize=12)
        # Hide the axes
        ax1.set_xticks([])
        ax1.set_yticks([])
        # Show the plot
        plt.tight_layout()
        fig1.savefig(f"DRAGON_RATES_{name_case}/{name_compo}/assembly_map_{fig_name}_g1.png", dpi=300)
        plt.close(fig1)

        ax2.set_xlim(0, grid_size)
        ax2.set_ylim(0, grid_size)
        ax2.set_aspect('equal')
        # Set the title and labels
        ax2.set_title(f"Assembly Map for {name_case} - {name_compo},thermal group", fontsize=16)
        ax2.set_xlabel("X-axis", fontsize=12)
        ax2.set_ylabel("Y-axis", fontsize=12)
        # Hide the axes
        ax2.set_xticks([])
        ax2.set_yticks([])
        # Show the plot
        plt.tight_layout()
        fig2.savefig(f"DRAGON_RATES_{name_case}/{name_compo}/assembly_map_{fig_name}_g2.png", dpi=300)
        plt.close(fig2)
        return


def plot_errors_BWR_assembly(errors_rates, name_case, name_compo, fig_name, unfold_symmetry, cmap="Spectral"):
    """
        Assembly map for ATRIUM-10

        C1 C2 C3 C5 C6 C5 C4 C3 C2 C1
        C2 C4 C7 C6 C7 C6 C6 C7 C4 C2
        C3 C7 C6 C6 C6 C7 C6 C6 C7 C3
        C5 C6 C6 C6 C6 C6 C7 C6 C5 C4
        C6 C7 C6 C6 WH WH WH C4 C6 C4
        C5 C6 C7 C6 WH WH WH C3 C7 C4
        C4 C6 C6 C7 WH WH WH C4 C4 C4
        C3 C7 C6 C6 C4 C3 C4 C4 C8 C3
        C2 C4 C7 C5 C6 C7 C4 C8 C4 C2
        C1 C2 C3 C4 C4 C4 C4 C3 C2 C1

        WH = Water Hole
        Cn = Fuel Cell number n

        Parameters:
        - errors_rates (list) : list containing dictionnaries with relative errors. Element in list corresponds to the group number.
            for each group number is a dictionary with key mix name and value the relative error.
        - name_case (str): Name of the geometry.
        - name_compo (str): Name of the MULTICOMPO obect with calculation parameter identifiers.
        - fig_name (str): Name of the figure to be saved.
        - unfold_symmetry (bool): If True, apply symmetry to the assembly map.
    
    """

    if unfold_symmetry:
            # Define the assembly map
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["C2", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["C3", "C7", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["C5", "C6", "C6", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["C6", "C7", "C6", "C6", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["C5", "C6", "C7", "C6", "WH", "WH", "WH", "C3", "C7", "C4"],
                ["C4", "C6", "C6", "C7", "WH", "WH", "WH", "C4", "C4", "C4"],
                ["C3", "C7", "C6", "C6", "C4", "C3", "C4", "C4", "C8", "C3"],
                ["C2", "C4", "C7", "C5", "C6", "C7", "C4", "C8", "C4", "C2"],
                ["C1", "C2", "C3", "C4", "C4", "C4", "C4", "C3", "C2", "C1"]
            ])
    else:
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["sy", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["sy", "sy", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["sy", "sy", "sy", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["sy", "sy", "sy", "sy", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "WH", "WH", "C3", "C7", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "WH", "C4", "C4", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C8", "C3"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C2"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C1"]
                ])

    a = os.path.exists(f"DRAGON_RATES_{name_case}/{name_compo}")
    print(f"DRAGON_RATES_{name_case}/{name_compo} already exists: {a}")
    if not a:
        os.makedirs(f"DRAGON_RATES_{name_case}/{name_compo}")
    for group_num in range(len(errors_rates)):
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(10, 10))

        # Define the grid size
        grid_size = 10

        min_error = min(list(errors_rates[group_num].values()))
        max_error = max(list(errors_rates[group_num].values()))
        # Normalize colors based on error values
        norm = plt.Normalize(min_error, max_error)
        color_map = plt.cm.get_cmap(cmap)

        # Loop through the assembly map and plot each cell
        for row in range(assembly_map.shape[0]):
            for col in range(assembly_map.shape[1]):
                cell_name = assembly_map[row][col]
                print(f"cell_name = {cell_name}")
                
                # Get the error on reaction rate for the current cell
                if cell_name in errors_rates[group_num].keys():
                    error_value = errors_rates[group_num][cell_name] if cell_name in errors_rates[group_num].keys() else 0
                    color = color_map(norm(error_value))  # Get color from colormap
                else:
                    rate_g1 = cell_name
                # Get the reaction rate for the current cell
                # Plot the cell with the corresponding color
                if cell_name == "sy" or cell_name == "WH":
                    print("skip") 
                else:
                    rect = plt.Rectangle((col, grid_size - row - 1), 1, 1, facecolor=color, edgecolor='black', linewidth=2)
                    ax.add_patch(rect)
                    # Add the reaction rate text inside the cell
                    ax.text(col + 0.5, grid_size - row - 0.5, f"{error_value:.2f}", ha='center', va='center', fontsize=12, color='black')
        # Set the limits and aspect ratio
        ax.set_xlim(0, grid_size)
        ax.set_ylim(0, grid_size)
        ax.set_aspect('equal')
        # Set the title and labels
        ax.set_title(f'{fig_name.replace("_"," ")}, group {group_num+1}', fontsize=16)
        # Hide the axes
        ax.set_xticks([])
        ax.set_yticks([])
        # Show the plot
        plt.tight_layout()
        fig.savefig(f"DRAGON_RATES_{name_case}/{name_compo}/assembly_map_{fig_name}_g{group_num+1}.png", dpi=300)
        plt.close(fig)

    return


def plot_pinwise_errors_BWR_assembly(errors_rates, name_case, name_compo, calculation_opt, fig_name, evaluation, cmap="bwr"):
    """
        Assembly map for ATRIUM-10

        C1 C2 C3 C5 C6 C5 C4 C3 C2 C1
        C2 C4 C7 C6 C7 C6 C6 C7 C4 C2
        C3 C7 C6 C6 C6 C7 C6 C6 C7 C3
        C5 C6 C6 C6 C6 C6 C7 C6 C5 C4
        C6 C7 C6 C6 WH WH WH C4 C6 C4
        C5 C6 C7 C6 WH WH WH C3 C7 C4
        C4 C6 C6 C7 WH WH WH C4 C4 C4
        C3 C7 C6 C6 C4 C3 C4 C4 C8 C3
        C2 C4 C7 C5 C6 C7 C4 C8 C4 C2
        C1 C2 C3 C4 C4 C4 C4 C3 C2 C1

        WH = Water Hole
        Cn = Fuel Cell number n

        Parameters:
        - errors_rates np.array with shape (num_groups, num_cells): array containing relative errors. Each row corresponds to a group number and each column to a cell.
        - name_case (str): Name of the geometry.
        - name_compo (str): Name of the MULTICOMPO obect with calculation parameter identifiers.
        - fig_name (str): Name of the figure to be saved.
        - unfold_symmetry (bool): If True, apply symmetry to the assembly map.
    
    """

    assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["sy", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["sy", "sy", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["sy", "sy", "sy", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["sy", "sy", "sy", "sy", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "WH", "WH", "C3", "C7", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "WH", "C4", "C4", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C8", "C3"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C2"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C1"]
                ])
    assembly_map = np.array([
                [10, 19, 27, 34, 40, 45, 49, 52, 54, 55],
                [ 9, 18, 26, 33, 39, 44, 48, 51, 53, 54],
                [ 8, 17, 25, 32, 38, 43, 47, 50, 51, 52],
                [ 7, 16, 24, 31, 00, 00, 00, 47, 48, 49],
                [ 6, 15, 23, 30, 00, 00, 00, 43, 44, 45],
                [ 5, 14, 22, 29, 00, 00, 00, 38, 39, 40],
                [ 4, 13, 21, 28, 29, 30, 31, 32, 33, 34],
                [ 3, 12, 20, 21, 22, 23, 24, 25, 26, 27],
                [ 2, 11, 12, 13, 14, 15, 16, 17, 18, 19],
                [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10],

                
                ])

    a = os.path.exists(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{name_compo}")
    print(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{name_compo} already exists: {a}")
    if not a:
        os.makedirs(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{name_compo}")

    ngroups = errors_rates.shape[0]
    ncells = errors_rates.shape[1]
    for gr in range(ngroups):
        
        error_grid = np.full_like(assembly_map, np.nan, dtype=float)

        for i in range(assembly_map.shape[0]):
            for j in range(assembly_map.shape[1]):
                pos = assembly_map[i, j]
                if pos != 0:
                    error_grid[i, j] = errors_rates[gr][pos-1]
        plt.figure(figsize=(10,10))
        im = plt.imshow(np.flipud(error_grid), origin="lower", cmap="bwr", vmin=np.min(errors_rates[gr]), vmax=np.max(errors_rates[gr]))

        plt.colorbar(im, label="Relative Difference (%)")
        
        nrows, ncols = error_grid.shape
        for i in range(nrows):
            for j in range(ncols):
                val = error_grid[i, j]
                if not np.isnan(val):
                    plt.text(j, nrows-1-i,f"{val:.2f}%",ha="center", va="center",color="black",fontsize=12,fontweight="bold")
        # Draw grid lines only inside the lattice
        ax = plt.gca()
        ax.set_xticks(np.arange(-0.5, ncols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, nrows, 1), minor=True)
        ax.grid(which="minor", color="black", linestyle="-", linewidth=0.8)

        # Hide the major ticks & labels completely
        ax.tick_params(which="both", bottom=False, left=False, labelbottom=False, labelleft=False)

        # Clip the grid strictly to the data area
        ax.set_xlim(-0.5, ncols-0.5)
        ax.set_ylim(-0.5, nrows-0.5)
                    
        if gr == 0:
            group_id = "thermal"
        elif gr == 1:
            group_id = "fast"
        if fig_name == "Fiss_over_abs_diff":
            plt.title(f'(D5-S2) Relative differences on $\\tau_f/\\tau_a$, {group_id} group', fontsize=16)
        else:
            plt.title(f'(D5-S2) Relative differences on fission rates, {group_id} group', fontsize=16)
        # Show the plot
        plt.tight_layout()
        plt.savefig(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{name_compo}/assembly_map_{fig_name}_g{gr+1}.png", dpi=300)
        plt.close()
        

    return


def plot_error_grid_from_list(plot_title, name_case, name_compo, group, error_list, N, cmap="bwr", text_color="black"):
    """
    Plots an N x N grid where each cell is colored according to the error value.
    
    Parameters:
    - name_case (str): Name of the geometry.
    - group (int): Group number.
    - error_list (list or numpy.ndarray): 1D list of error values, ordered by increasing x first, then increasing y.
    - N (int): Size of the grid (assumed to be N x N).
    - cmap (str): Colormap name (default is 'coolwarm').
    - text_color (str): Color of the text inside the squares (default is black).
    
    Returns:
    - None (displays the plot)
    """
    # Convert list into an NxN array
    error_matrix = np.array(error_list).reshape(N, N)  # Reshape row-wise

    fig, ax = plt.subplots(figsize=(8, 8))

    # Normalize colors based on error values
    min_err, max_err = np.min(error_matrix), np.max(error_matrix)
    norm = plt.Normalize(min_err, max_err)
    color_map = plt.cm.get_cmap(cmap)

    # Plot each square
    for i in range(N):
        for j in range(N):
            value = error_matrix[i, j]
            color = color_map(norm(value))  # Get color from colormap
            
            # The input ordering assumes row-major format (increasing x first), 
            # but Matplotlib's (0,0) is at the top-left, so we flip vertically
            ax.add_patch(plt.Rectangle((j, i), 1, 1, color=color, edgecolor='black'))

            # Add text at the center of each square
            ax.text(j + 0.5, i + 0.5, f"{value:.2f} %", 
                    ha='center', va='center', fontsize=15, color=text_color)

    # Set axis limits
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)

    # Remove axes ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=color_map, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Relative errors on $\\tau_f$ for U235 (%)", fontsize=12)

    # Show the plot
    plt.title("Spatial Error (D5-S2) Distribution on U235 fission rates", fontsize=14)
    plt.savefig(f"DRAGON_RATES_{name_case}/{plot_title}.png", dpi=300)
    plt.close(fig)
    return


def plot_spectrum_comparison(energy_mesh, FLUX_groups, S2_spectra, name_case, CPO_name, calculation_opt, evaluation):
    n_groups = len(FLUX_groups)
    try:
        S2_spectrum = S2_spectra[f"{n_groups}g"]
    except KeyError:
        print(f"Error: Serpent2 spectrum for {n_groups} groups not found in S2_spectra dictionary.")
        return    
    
    a = os.path.exists(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}")
    print(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name} already exists: {a}")
    if not a:
        os.makedirs(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}")
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
    plt.savefig(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}/flux_spectrum_comparison.png")
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
    plt.savefig(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}/flux_spectrum_comparison_lethargy.png")
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
    plt.savefig(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}/relative_difference_flux_spectrum.png")
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
    plt.savefig(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}/absolute_difference_flux_spectrum.png")
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
    plt.savefig(f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}/difference_flux_spectrum.png")
    plt.close()

    
def plot_U238_abs_rates(energy_mesh, D5_U238_abs_rates, S2_U238_abs_rates, name_case, CPO_name, calculation_opt, evaluation):
    """ 
    plot spatially integrated U238 absorption rates on energy mesh
    """
    results_dir = f"DRAGON_RATES_{name_case}_{evaluation}/{calculation_opt}/{CPO_name}"
    a = os.path.exists(results_dir)
    print(f"{results_dir} already exists: {a}")
    if not a:
        os.makedirs(results_dir)
        
    # Normalize DRAGON rates to Serpent2 rates for better comparison
    normalization_factor = np.sum(S2_U238_abs_rates) / np.sum(D5_U238_abs_rates)
    D5_U238_abs_rates = D5_U238_abs_rates * normalization_factor
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