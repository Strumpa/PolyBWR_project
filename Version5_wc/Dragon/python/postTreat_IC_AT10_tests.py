import re
from typing import Dict, List, Tuple
from postTreat_IC_anis_tests import *
 
if __name__ == "__main__":
    series_of_tests = "AT10_GLOW_tests"
    # Example usage
    result_file = f"../Linux_aarch64/{series_of_tests}.result"
    
    try:
        tests = parse_ic_anis_tests(result_file)
        print_test_summary(tests)
        
        # Export to CSV
        print("\n")
        export_to_csv(tests, f"{series_of_tests}_results.csv")
        
        # Analyze keff differences (in pcm)
        analyze_keff_differences(tests)
        
        # Compare methods for GEO geometry
        test_names = ["AT10_2x2", 
                      "AT10_3x3", "AT10_3x3Gd",
                      "AT10_4x4", "AT10_4x4Gd",
                      "AT10_LAT", "AT10_LAT_Gd",
                      "AT10_LAT_MODEBOX", "AT10_LAT_MODEBOX_Gd",
                      "AT10_A", "AT10_A_Gd"]
        geometry_generators = ['GEO', 'GLOW']
        for test_name in test_names:
            # Compare GEO vs GLOW for all solvers
            print("\n")
            compare_geo_vs_glow(tests, test_name)
            for geometry in geometry_generators:
                # Compare methods for GLOW geometry
                print("\n")
                compare_methods_by_geometry(tests, test_name, geometry=geometry)
        
        # Compare GEO vs GLOW for all tests
        print("\n")
        compare_geo_vs_glow(tests)
        
    except FileNotFoundError:
        print(f"Error: File '{result_file}' not found.")
    except Exception as e:
        print(f"Error parsing file: {e}")

