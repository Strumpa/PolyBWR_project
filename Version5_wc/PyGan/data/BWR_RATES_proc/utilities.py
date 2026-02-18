# Utility functions for processing BWR rates data
from collections import defaultdict

def sum_rates_over_iso(rates_dict):
    """
    rates_dict is a nested dictionary with the structure:
    {
        "iso1": {
            "mix1": [val1, val2, ...],
            "mix2": [...],
            ...
        },
        ...
    }

    This function sums the values over isotopes while preserving the mix structure,
    returning:
    {
        "mix1": [sum_val1, sum_val2, ...],
        "mix2": [...],
        ...
    }
    """
    result = defaultdict(list)

    for isotope_data in rates_dict.values():
        for mix, values in isotope_data.items():
            if mix not in result:
                result[mix] = [0.0] * len(values)
            for i, val in enumerate(values):
                result[mix][i] += val

    return dict(result)