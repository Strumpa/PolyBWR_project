import re

def extract_numbers_from_list(input_list):
    # Initialize lists to store extracted numbers, materials, densities, and additional info
    numbers = []
    materials = []
    densities = []
    additional_info = []

    # Define regex pattern to match numbers, materials, densities, and additional info
    pattern = r'(?:#?\d+)\s+([-+]?\d+(\.\d+)?)\s+([-+]?\d+(\.\d+)?)\s+imp:n=1\s+\$(.*)'

    for item in input_list:
        # Find all matches in the string using regex
        matches = re.match(pattern, item)
        
        if matches:
            # Extract numbers, materials, densities, and additional info
            number = float(matches.group(1))
            material = int(matches.group(2))
            density = float(matches.group(3))
            info = matches.group(4).strip()  # Extract and strip whitespace from additional info

            # Append extracted values to respective lists
            numbers.append(number)
            materials.append(material)
            densities.append(density)
            additional_info.append(info)
        else:
            print(f"No match found for item: {item}")

    return numbers, materials, densities, additional_info

# Example usage:
input_list = ['#25 23 -10.461000  -40  38 -39   imp:n=1  $ FUE m23',
              '2 22  -0.001000  40 -41  38 -39   imp:n=1  $ HEL m22',
              '3  5  -6.550000  41 -42  38 -39   imp:n=1  $ CAN m05',
              '4  2  -0.458426  5  -6  22 -23  42  38 -39   imp:n=1  $ COO m02']

numbers, materials, densities, additional_info = extract_numbers_from_list(input_list)

# Print the extracted numbers, materials, densities, and additional info
print("Numbers:", numbers)
print("Materials:", materials)
print("Densities:", densities)
print("Additional Info:", additional_info)
