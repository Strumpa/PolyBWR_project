import re
from typing import Dict, List, Tuple

def parse_ic_anis_tests(file_path: str) -> List[Dict]:
    """
    Parse the IC_anis_tests.result file and extract test information.
    
    For each test, extracts:
    - Test name and keff
    - SALT, ASM and FLU timing information for different solvers
    
    Args:
        file_path: Path to the IC_anis_tests.result file
        
    Returns:
        List of dictionaries containing test information
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    tests = []
    
    # Find all test completion messages
    # Pattern: >|test multicell {name} {method} completed, keff=  {keff}
    test_pattern = r'>+\|test multicell ([^ ]+) ([^,]+) completed,\s*keff=\s*([e\d\.\-\+]+)'
    test_matches = list(re.finditer(test_pattern, content))
    
    # Track SALT tracking information by test name (geometry)
    # All 4 methods for a geometry share the same 3 SALT tracking invocations
    # SALT[0] = IC tracking (for PIJ/FLUX-CURRENT)
    # SALT[1] = CP tracking (for CP)
    # SALT[2] = MOC tracking (for MOC)
    salt_tracking_cache = {}
    
    for match in test_matches:
        test_name = match.group(1)
        method = match.group(2).strip()
        keff = float(match.group(3))
        test_start_pos = match.start()
        
        # Find the next test or end of file
        if test_matches.index(match) + 1 < len(test_matches):
            test_end_pos = test_matches[test_matches.index(match) + 1].start()
        else:
            test_end_pos = len(content)
        
        # Extract the section for this test (look backwards from test completion)
        # We need to find the ASM, FLU, and SALT sections that precede this test
        section_start = max(0, test_start_pos - 100000)  # Look back 100k characters
        section_content = content[section_start:test_start_pos]
        
        # Find ASM and FLU module information (most recent)
        asm_info = extract_module_info(section_content, 'ASM', occurrence=-1)
        flu_info = extract_module_info(section_content, 'FLU', occurrence=-1)
        
        # For SALT tracking:
        # All 3 SALT trackings (IC, CP, MOC) happen before the FULL PIJ test completes
        # They are then reused by FLUX-CURRENT, CP, and MOC tests respectively
        # We cache all 3 when processing FULL PIJ, then retrieve them for other methods
        
        if 'FULL PIJ' in method:
            # Extract all 3 SALT trackings from the section
            salt_ic = extract_module_info(section_content, 'SALT', occurrence=0)
            salt_cp = extract_module_info(section_content, 'SALT', occurrence=1)
            salt_moc = extract_module_info(section_content, 'SALT', occurrence=-1)
            
            # Cache all 3 for this test geometry
            salt_tracking_cache[test_name] = {
                'IC': salt_ic,
                'CP': salt_cp,
                'MOC': salt_moc
            }
            
            # FULL PIJ uses IC tracking
            salt_info = salt_ic
            
        elif 'FLUX-CURRENT' in method:
            # Reuse IC tracking from cache
            if test_name in salt_tracking_cache:
                salt_info = salt_tracking_cache[test_name]['IC']
            else:
                salt_info = {}
                
        elif 'CP' in method:
            # Reuse CP tracking from cache
            if test_name in salt_tracking_cache:
                salt_info = salt_tracking_cache[test_name]['CP']
            else:
                salt_info = {}
                
        elif 'MOC' in method:
            # Reuse MOC tracking from cache
            if test_name in salt_tracking_cache:
                salt_info = salt_tracking_cache[test_name]['MOC']
            else:
                salt_info = {}
        else:
            # Fallback
            salt_info = {}
        
        test_entry = {
            'name': test_name,
            'method': method,
            'keff': keff,
            'salt_time': salt_info.get('time'),
            'salt_memory': salt_info.get('memory'),
            'asm_time': asm_info.get('time'),
            'asm_memory': asm_info.get('memory'),
            'flu_time': flu_info.get('time'),
            'flu_memory': flu_info.get('memory'),
            'flu_kinf': flu_info.get('kinf'),  # Final KINF from flux calculation
        }
        
        tests.append(test_entry)
    
    return tests


def extract_module_info(content: str, module_name: str, occurrence: int = -1) -> Dict:
    """
    Extract timing and memory information for a specific module from the content.
    
    Args:
        content: Text content to search
        module_name: Name of the module (e.g., 'ASM', 'FLU', 'SALT')
        occurrence: Which occurrence to extract (0=first, -1=last, -2=second-to-last, -3=third-to-last)
        
    Returns:
        Dictionary with 'time', 'memory', and optionally 'kinf' keys
    """
    info = {}
    
    # Pattern for module timing: -->>MODULE {name}:        : TIME SPENT=        {time} MEMORY USAGE= {memory}
    # The name may be followed by spaces, then a colon and more spaces
    timing_pattern = rf'-->>MODULE\s+{module_name}:\s*:\s*TIME SPENT=\s*([\d\.E\+\-]+)\s+MEMORY USAGE=\s*([\d\.E\+\-]+)'
    timing_match = list(re.finditer(timing_pattern, content))
    
    if timing_match:
        try:
            # Get the specified occurrence
            match = timing_match[occurrence]
            info['time'] = float(match.group(1))
            info['memory'] = float(match.group(2))
        except IndexError:
            # Occurrence not found, leave info empty
            pass
    
    # For FLU module, also extract KINF value
    if module_name == 'FLU':
        # Pattern: FINAL KINF= {kinf}
        kinf_pattern = r'FINAL KINF=\s*([\d\.E\+\-]+)'
        kinf_match = list(re.finditer(kinf_pattern, content))
        if kinf_match:
            info['kinf'] = float(kinf_match[-1].group(1))
    
    return info


def print_test_summary(tests: List[Dict]) -> None:
    """
    Print a formatted summary of all parsed tests.
    
    Args:
        tests: List of test dictionaries from parse_ic_anis_tests()
    """
    print("\n" + "="*140)
    print("IC ANIS TESTS SUMMARY")
    print("="*140)
    print(f"{'Test Name':<25} {'Method':<28} {'Keff':<15} {'SALT':<10} {'ASM':<10} {'FLU':<10} {'FLU KINF':<15}")
    print("-"*140)
    
    for test in tests:
        test_name = test['name']
        method = test['method'][:26] if test['method'] else "N/A"
        keff = f"{test['keff']:.6e}" if test['keff'] else "N/A"
        salt_time = f"{test['salt_time']:.3f}" if test['salt_time'] is not None else "N/A"
        asm_time = f"{test['asm_time']:.3f}" if test['asm_time'] is not None else "N/A"
        flu_time = f"{test['flu_time']:.3f}" if test['flu_time'] is not None else "N/A"
        flu_kinf = f"{test['flu_kinf']:.6e}" if test['flu_kinf'] is not None else "N/A"
        
        print(f"{test_name:<25} {method:<28} {keff:<15} {salt_time:<10} {asm_time:<10} {flu_time:<10} {flu_kinf:<15}")
    
    print("="*140 + "\n")


def get_tests_by_name(tests: List[Dict], test_name: str) -> List[Dict]:
    """
    Filter tests by name.
    
    Args:
        tests: List of test dictionaries
        test_name: Name of test to filter
        
    Returns:
        Filtered list of tests
    """
    return [t for t in tests if t['name'] == test_name]


def get_tests_by_method(tests: List[Dict], method: str) -> List[Dict]:
    """
    Filter tests by method.
    
    Args:
        tests: List of test dictionaries
        method: Method name to filter (e.g., 'FULL PIJ reconstruction')
        
    Returns:
        Filtered list of tests
    """
    return [t for t in tests if method.lower() in t['method'].lower()]


def get_geometry_type(test: Dict) -> str:
    """
    Determine the geometry generator type (GEO or GLOW).
    
    Args:
        test: Test dictionary
        
    Returns:
        'GLOW' if method contains '(GLOW)', otherwise 'GEO'
    """
    return 'GLOW' if '(GLOW)' in test['method'] or 'GLOW' in test['method'] else 'GEO'


def get_base_method(test: Dict) -> str:
    """
    Extract the base method name without geometry generator prefix.
    
    Args:
        test: Test dictionary
        
    Returns:
        Base method name (e.g., 'FULL PIJ reconstruction' from '(GLOW) FULL PIJ reconstruction')
    """
    method = test['method']
    # Remove (GLOW) prefix if present
    if method.startswith('(GLOW) '):
        return method[7:].strip()
    if method.startswith('(GLOW+MACRO) '):
        return method[12:].strip()
    return method.strip()


def compare_methods_by_geometry(tests: List[Dict], test_name: str, geometry: str = 'GEO', reference_method: str = 'MOC') -> None:
    """
    Compare different methods for a specific test and geometry generator.
    Differences are expressed in pcm relative to a reference method.
    
    Args:
        tests: List of test dictionaries
        test_name: Name of the test to compare
        geometry: Geometry generator ('GEO' or 'GLOW')
        reference_method: Reference method for comparison (default: 'FULL PIJ reconstruction')
    """
    # Filter tests by name and geometry
    test_data = [t for t in tests if t['name'] == test_name and get_geometry_type(t) == geometry]
    
    if not test_data:
        print(f"No tests found for {test_name} with {geometry} geometry")
        return
    
    # Find reference test
    ref_test = None
    for t in test_data:
        if get_base_method(t) == reference_method:
            ref_test = t
            break
    
    if not ref_test:
        print(f"Warning: Reference method '{reference_method}' not found. Using first test as reference.")
        ref_test = test_data[0]
    
    ref_keff = ref_test['keff']
    
    print(f"\n{'='*115}")
    print(f"METHOD COMPARISON FOR {test_name} - {geometry} GEOMETRY")
    print(f"Reference: {get_base_method(ref_test)} (keff = {ref_keff:.6e})")
    print(f"{'='*115}")
    print(f"{'Method':<35} {'Keff':<15} {'Diff (pcm)':<15} {'SALT':<10} {'ASM':<10} {'FLU':<10}")
    print(f"{'-'*115}")
    
    for test in test_data:
        base_method = get_base_method(test)
        keff = test['keff']
        diff_pcm = (keff - ref_keff) * 1e5
        salt_time = f"{test['salt_time']:.3f}" if test['salt_time'] is not None else "N/A"
        asm_time = f"{test['asm_time']:.3f}" if test['asm_time'] is not None else "N/A"
        flu_time = f"{test['flu_time']:.3f}" if test['flu_time'] is not None else "N/A"
        
        marker = " (REF)" if test == ref_test else ""
        print(f"{base_method:<35} {keff:<15.6e} {diff_pcm:<15.2f} {salt_time:<10} {asm_time:<10} {flu_time:<10}{marker}")


def compare_geo_vs_glow(tests: List[Dict], test_name: str = None) -> None:
    """
    Compare GEO vs GLOW results for each solver method.
    Differences are expressed in pcm: (keff_GLOW - keff_GEO) * 1e5
    
    Args:
        tests: List of test dictionaries
        test_name: Optional - compare for specific test only. If None, compare all tests.
    """
    # Group tests by name and base method
    if test_name:
        test_names = [test_name]
        test_subset = get_tests_by_name(tests, test_name)
    else:
        test_names = sorted(set(t['name'] for t in tests))
        test_subset = tests
    
    print(f"\n{'='*120}")
    print("GEO vs GLOW COMPARISON (Diff in pcm: keff_GLOW - keff_GEO)")
    print(f"{'='*120}")
    print(f"{'Test':<15} {'Method':<30} {'GEO Keff':<15} {'GLOW Keff':<15} {'Diff (pcm)':<15}")
    print(f"{'-'*120}")
    
    for name in test_names:
        name_tests = [t for t in test_subset if t['name'] == name]
        
        # Group by base method
        methods = {}
        for test in name_tests:
            base_method = get_base_method(test)
            geo_type = get_geometry_type(test)
            
            if base_method not in methods:
                methods[base_method] = {}
            methods[base_method][geo_type] = test
        
        # Compare GEO vs GLOW for each method
        for method in sorted(methods.keys()):
            if 'GEO' in methods[method] and 'GLOW' in methods[method]:
                geo_test = methods[method]['GEO']
                glow_test = methods[method]['GLOW']
                
                geo_keff = geo_test['keff']
                glow_keff = glow_test['keff']
                diff_pcm = (glow_keff - geo_keff) * 1e5
                
                print(f"{name:<15} {method:<30} {geo_keff:<15.6e} {glow_keff:<15.6e} {diff_pcm:<15.2f}")


def export_to_csv(tests: List[Dict], output_file: str) -> None:
    """
    Export parsed tests to a CSV file.
    
    Args:
        tests: List of test dictionaries
        output_file: Path to output CSV file
    """
    import csv
    
    try:
        with open(output_file, 'w', newline='') as f:
            fieldnames = ['name', 'method', 'keff', 'salt_time', 'salt_memory', 'asm_time', 'asm_memory', 'flu_time', 'flu_memory', 'flu_kinf']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            
            writer.writeheader()
            for test in tests:
                writer.writerow(test)
        
        print(f"Exported {len(tests)} tests to {output_file}")
    except Exception as e:
        print(f"Error exporting to CSV: {e}")


def compare_methods_for_test(tests: List[Dict], test_name: str) -> None:
    """
    Compare results for different methods applied to the same test.
    
    Args:
        tests: List of test dictionaries
        test_name: Name of the test to compare
    """
    test_data = get_tests_by_name(tests, test_name)
    
    if not test_data:
        print(f"No tests found with name '{test_name}'")
        return
    
    print(f"\n{'='*100}")
    print(f"COMPARISON FOR TEST: {test_name}")
    print(f"{'='*100}")
    print(f"{'Method':<35} {'Keff':<15} {'ASM Time':<12} {'FLU Time':<12} {'FLU KINF':<15}")
    print(f"{'-'*100}")
    
    for test in test_data:
        method = test['method'][:33] if test['method'] else "N/A"
        keff = f"{test['keff']:.6e}" if test['keff'] else "N/A"
        asm_time = f"{test['asm_time']:.3f}" if test['asm_time'] is not None else "N/A"
        flu_time = f"{test['flu_time']:.3f}" if test['flu_time'] is not None else "N/A"
        flu_kinf = f"{test['flu_kinf']:.6e}" if test['flu_kinf'] is not None else "N/A"
        
        print(f"{method:<35} {keff:<15} {asm_time:<12} {flu_time:<12} {flu_kinf:<15}")


def analyze_keff_differences(tests: List[Dict]) -> None:
    """
    Analyze differences in keff between different methods for each test.
    Differences are expressed in pcm (percent mille): (keff_test - keff_ref) * 1e5
    
    Args:
        tests: List of test dictionaries
    """
    # Group tests by name
    test_names = set(t['name'] for t in tests)
    
    print(f"\n{'='*100}")
    print("KEFF DIFFERENCES ANALYSIS BY TEST (in pcm)")
    print(f"{'='*100}")
    print(f"{'Test Name':<20} {'Min Keff':<15} {'Max Keff':<15} {'Diff (pcm)':<15}")
    print(f"{'-'*100}")
    
    for name in sorted(test_names):
        test_data = get_tests_by_name(tests, name)
        keffs = [t['keff'] for t in test_data]
        
        if len(keffs) > 1:
            min_keff = min(keffs)
            max_keff = max(keffs)
            diff_pcm = (max_keff - min_keff) * 1e5
            
            print(f"{name:<20} {min_keff:<15.6e} {max_keff:<15.6e} {diff_pcm:<15.2f}")


def compare_test_pairs(tests: List[Dict], test_pairs: List[Tuple[str, str]]) -> None:
    """
    Compare keff values between pairs of tests for all geometry and method combinations.
    Differences are expressed in pcm (percent mille): (keff_test2 - keff_test1) * 1e5
    
    This function allows checking that certain test pairs give the same results
    across all combinations of geometry generators (GEO, GLOW) and methods 
    (FULL PIJ, FLUX-CURRENT, CP, MOC).
    
    Args:
        tests: List of test dictionaries from parse_ic_anis_tests()
        test_pairs: List of tuples containing test name pairs to compare
                    e.g., [('3x3_Gd_BL', '3x3_Gd_TR'), ('3x3_Gd_BC', '3x3_Gd_TC')]
    
    Example:
        compare_test_pairs(tests, [
            ('3x3_Gd_BL', '3x3_Gd_TR'),  # Bottom-Left vs Top-Right
            ('3x3_Gd_BC', '3x3_Gd_TC')   # Bottom-Center vs Top-Center
        ])
    """
    # Define all possible combinations
    geometries = ['GEO', 'GLOW']
    methods = ['FULL PIJ reconstruction', 'FLUX-CURRENT iteration', 'CP', 'MOC']
    
    for test1_name, test2_name in test_pairs:
        print(f"\n{'='*120}")
        print(f"COMPARISON: {test1_name} vs {test2_name}")
        print(f"Difference shown as: (keff_{test2_name} - keff_{test1_name}) × 10⁵ [pcm]")
        print(f"{'='*120}")
        print(f"{'Geometry':<12} {'Method':<30} {f'{test1_name} keff':<18} {f'{test2_name} keff':<18} {'Diff (pcm)':<15} {'Status':<10}")
        print(f"{'-'*120}")
        
        # Track statistics
        all_diffs = []
        
        for geometry in geometries:
            for method in methods:
                # Get the tests for both test names
                test1_data = [t for t in tests if t['name'] == test1_name and 
                             get_geometry_type(t) == geometry and 
                             get_base_method(t) == method]
                test2_data = [t for t in tests if t['name'] == test2_name and 
                             get_geometry_type(t) == geometry and 
                             get_base_method(t) == method]
                
                if test1_data and test2_data:
                    test1 = test1_data[0]
                    test2 = test2_data[0]
                    
                    keff1 = test1['keff']
                    keff2 = test2['keff']
                    diff_pcm = (keff2 - keff1) * 1e5
                    all_diffs.append(abs(diff_pcm))
                    
                    # Determine status (flagging large differences)
                    if abs(diff_pcm) < 0.1:
                        status = "✓"
                    elif abs(diff_pcm) < 1.0:
                        status = "~"
                    else:
                        status = "⚠"
                    
                    # Shorten method name for display
                    method_short = method[:28] if len(method) > 28 else method
                    
                    print(f"{geometry:<12} {method_short:<30} {keff1:<18.6e} {keff2:<18.6e} {diff_pcm:<15.2f} {status:<10}")
                else:
                    print(f"{geometry:<12} {method:<30} {'N/A':<18} {'N/A':<18} {'N/A':<15} {'?':<10}")
        
        # Print summary statistics
        if all_diffs:
            print(f"{'-'*120}")
            print(f"Statistics: Max |diff| = {max(all_diffs):.2f} pcm, "
                  f"Mean |diff| = {sum(all_diffs)/len(all_diffs):.2f} pcm, "
                  f"Min |diff| = {min(all_diffs):.2f} pcm")
            print(f"Legend: ✓ = |diff| < 0.1 pcm (excellent), ~ = |diff| < 1 pcm (good), ⚠ = |diff| ≥ 1 pcm (review)")


if __name__ == "__main__":
    series_of_tests = "IC_anis_MACRO_tests"
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
        test_names = ["3x3_Gd_C", 
                      "3x3_Gd_TR", "3x3_Gd_BL",
                      "3x3_Gd_TC", "3x3_Gd_BC",
                      "3x3_Gd_LC", "3x3_Gd_RC",
                      "3x3_UOX_C",
                      "3x3_UOX_TR", "3x3_UOX_BL",
                      "3x3_UOX_TC", "3x3_UOX_BC"]
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
        
        # Compare test pairs across all geometry/method combinations
        print("\n")
        compare_test_pairs(tests, [
            ('3x3_Gd_BL', '3x3_Gd_TR'),  # Bottom-Left vs Top-Right (symmetry check)
            ('3x3_Gd_BC', '3x3_Gd_TC'),   # Bottom-Center vs Top-Center (symmetry check)
            ('3x3_Gd_RC', '3x3_Gd_LC'),    # Bottom-Right vs Top-Left (symmetry check)
            ('3x3_UOX_BL', '3x3_UOX_TR'),  # Bottom-Left vs Top-Right (symmetry check)
            ('3x3_UOX_BC', '3x3_UOX_TC'),   # Bottom-Center vs Top-Center (symmetry check)
        ])
        
    except FileNotFoundError:
        print(f"Error: File '{result_file}' not found.")
    except Exception as e:
        print(f"Error parsing file: {e}")

