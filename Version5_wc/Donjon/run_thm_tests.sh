#!/bin/bash

# bash script running a list of commands for THM module tests
# Usage: ./run_thm_tests.sh
# 2 Power normalizations : 38.4kW or 9.6kW, 2 axial discretizations : 10 or 20 control volumes, power distributions : flat or sine
# 2 options for pressure drop : 0 or 1, DFM not implemented.
# Create df for each test case to be used as non regresion of perssure drop model before implementation of DFM

# Define the list of test cases to run
test_cases=(
    "thm_10_flat_pow_9600W_PDROP0_DFM0"
    "thm_10_flat_pow_9600W_PDROP1_DFM0"
    "thm_10_flat_pow_38400W_PDROP0_DFM0"
    "thm_10_flat_pow_38400W_PDROP1_DFM0"

    "thm_20_flat_pow_9600W_PDROP0_DFM0"
    "thm_20_flat_pow_9600W_PDROP1_DFM0"
    "thm_20_flat_pow_38400W_PDROP0_DFM0"
    "thm_20_flat_pow_38400W_PDROP1_DFM0"

    "thm_10_sine_pow_9600W_PDROP0_DFM0"
    "thm_10_sine_pow_9600W_PDROP1_DFM0"
    "thm_10_sine_pow_38400W_PDROP0_DFM0"
    "thm_10_sine_pow_38400W_PDROP1_DFM0"

    "thm_20_sine_pow_9600W_PDROP0_DFM0"
    "thm_20_sine_pow_9600W_PDROP1_DFM0"
    "thm_20_sine_pow_38400W_PDROP0_DFM0"
    "thm_20_sine_pow_38400W_PDROP1_DFM0"
)

# Define the path to the post-processing script
python_post_treatment_path="./Linux_aarch64/parse_and_create_thm_df.py"

# Loop through each test case
for case in "${test_cases[@]}"; do
    echo "Running simulation: $case"
    ./rdonjon "$case.x2m"
    if [ $? -ne 0 ]; then
        echo "Error: Simulation failed for $case"
        exit 1
    fi

    echo "Post-processing: $case.result"
    python3 "$python_post_treatment_path" "./Linux_aarch64/$case.result"
    if [ $? -ne 0 ]; then
        echo "Error: Post-processing failed for $case"
        exit 1
    fi

    echo "Test completed successfully: $case"
    echo "----------------------------------------"
done

echo "All tests completed successfully."
