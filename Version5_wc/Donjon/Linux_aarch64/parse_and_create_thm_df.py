import pandas as pd
import re
import matplotlib.pyplot as plt
import os
import sys

def parse_results_table(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Filter out separator lines
    content_lines = [line for line in lines if not re.match(r'^\s*[|_]+\s*$', line)]

    # Extract header and unit lines
    header_line = content_lines[0]
    unit_line = content_lines[1]

    # Extract header names and units between '|', skip empty ends
    col_names = [col.strip() for col in header_line.strip().split('|')[1:-1]]
    col_units = [unit.strip() for unit in unit_line.strip().split('|')[1:-1]]

    # Combine names and units
    headers = [
        f"{name} ({unit})" if name and unit else name or unit or f"Col{i}"
        for i, (name, unit) in enumerate(zip(col_names, col_units))
    ]

    # Extract data rows
    data_rows = []
    for line in content_lines[2:]:
        fields = [f.strip() for f in line.strip().split('|')[1:-1]]
        if len(fields) == len(headers):
            data_rows.append(fields)
        else:
            print(f"⚠️ Skipping malformed line ({len(fields)} fields): {line.strip()}")

    # Build DataFrame
    df = pd.DataFrame(data_rows, columns=headers)

    # Convert columns to numeric where possible
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='ignore')  # keep 'TOP', 'BOT' etc. as strings

    return df

def parse_result_file(filename, n_nodes):
    """
    Parses a DONJON .result file and returns a DataFrame with the data.
    """
    results_df_name = f"df_test_{filename.replace('.result', '.csv')}"   
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Initialize an empty list to store the additional data to be stored in the DataFrame
    convective_heat_transfer_coeffs = []    
    T_surface_clad = []  # List to store the clad surface temperatures
    for line in lines:
        ## Add a condition if HCONV or TSCLAD is found
        if "HCONV =" in line: # retrieve convective heat transfer coef in W/m^2/K
            print(line.strip())
            # Extract the convective heat transfer coefficient
            hconv = float(line.strip().split('=')[1].strip())
            convective_heat_transfer_coeffs.append(hconv)
        elif "TSCLAD =" in line: # retrieve clad surface temperature in K
            print(line.strip())
            # Extract the clad surface temperature
            tsclad = float(line.strip().split('=')[1].strip())
            T_surface_clad.append(tsclad)
        if line.strip() == "|     |   TCOMB    |   TSURF    |    DCOOL    |    TCOOL    |    PCOOL    |    HCOOL    |    QFUEL    |    QCOOL    |    VOID   |     QUAL    |     SLIP    |  FLOW  |":
            print("Begin parsing THM table")
            #print(lines[lines.index(line)+3:lines.index(line)+3+n_nodes])
            table = lines[lines.index(line)+3:lines.index(line)+3+n_nodes]
            break
    # Store the table data in a list
    for i in range(len(table)):
        # Remove leading and trailing whitespace and split by '|'
        table[i] = table[i].strip().split('|')
        # Remove empty strings from the list
        table[i] = [x.strip() for x in table[i] if x.strip()]

        print(f"Row {i}: {table[i]}")

    # Create a DataFrame from the table
    df = pd.DataFrame(table, columns=['Z Mesh', 'TCOMB (K)', 'TSURF (K)', 'DCOOL (Kg/m3)', 'TCOOL (K)',
                                      'PCOOL (Pa)', 'HCOOL (J/Kg)', 'QFUEL (W/m3)', 'QCOOL (W/m3)',
                                      'VOID', 'QUAL', 'SLIP', 'FLOW'])
    
    print(f"DataFrame created with {len(df)} rows and {len(df.columns)} columns.")
    if len(convective_heat_transfer_coeffs) > 0:
        # Add the convective heat transfer coefficients to the DataFrame
        df['HCONV (W/m2/K)'] = convective_heat_transfer_coeffs[-n_nodes::]  # Use the last n_nodes values
    if len(T_surface_clad) > 0:
        # Add the clad surface temperatures to the DataFrame
        df['TSCLAD (K)'] = T_surface_clad[-n_nodes::]
    
    # save df to csv file in dataFrames directory
    if not os.path.exists(f'{os.environ["HOME"]}/working_dir/BWR-multiphysics/Version5_ev3785/Donjon/Linux_aarch64/dataFrames_wc'):
        os.makedirs(f'{os.environ["HOME"]}/working_dir/BWR-multiphysics/Version5_ev3785/Donjon/Linux_aarch64/dataFrames_wc')
        
    df.to_csv(f"{os.environ['HOME']}/working_dir/BWR-multiphysics/Version5_ev3785/Donjon/Linux_aarch64/dataFrames_wc/{results_df_name}", index=False)

    return


def parse_results_name(results_name):
    """
    Parses the results filename and extracts parameters based on the expected pattern:
    "thm_{nz}_h{height:.0f}_{shape}_pow_{Ptot}W_{correl}_PDROP{PDROP_opt}_DFM{DFM_opt}"
    
    Returns:
        dict: Dictionary containing extracted values.
    """
    pattern = (
        r"thm_(?P<nz>\d+)_h(?P<height>\d+)_"
        r"(?P<shape>\w+)_pow_(?P<Ptot>\d+)W_"
        r"(?P<correl>\w+)_PDROP(?P<PDROP_opt>\d+)_DFM(?P<DFM_opt>\d+)"
    )
    
    match = re.match(pattern, results_name)
    if not match:
        raise ValueError("Filename does not match expected format.")
    
    # Convert numeric values to int
    result = match.groupdict()
    result['nz'] = int(result['nz'])
    result['height'] = int(result['height'])
    result['Ptot'] = int(result['Ptot'])
    result['shape'] = result['shape']  # Keep as string
    result['correl'] = result['correl']  # Keep as string
    result['PDROP_opt'] = int(result['PDROP_opt'])
    result['DFM_opt'] = int(result['DFM_opt'])

    return result


if __name__ == "__main__":
    default_path_to_results  = "thm_10_sine_pow_38400W_PDROP0_DFM0.result"
    if len(sys.argv) > 1:
        pwd_path = os.getcwd()  # Get the current working directory
        results_path = sys.argv[1] 
        results_name = results_path.split("/")[-1]  # Get the last part of the path ie file name
        change_to_dir = "Linux_aarch64"  # Change to the directory where the results file is located
        os.chdir(change_to_dir)  # Change to the directory of the results file
    else: 
        results_name = default_path_to_results
        change_to_dir = None
    
    result_options = parse_results_name(results_name)
    # Extract parameters from the results name
    z_nodes = result_options['nz']  # Number of axial mesh points
    height = result_options['height']  # Height of the fuel assembly
    power_shape = result_options['shape']  # Shape of the fuel assembly
    total_power = result_options['Ptot']  # Total power in W
    correlation = result_options['correl']  # Power profile type (flat/sine)
    PDROP_opt = result_options['PDROP_opt']  # Pressure drop model option (0/1)
    DFM_opt = result_options['DFM_opt']  # DFM model option (0/1)
    parse_result_file(results_name, z_nodes)
