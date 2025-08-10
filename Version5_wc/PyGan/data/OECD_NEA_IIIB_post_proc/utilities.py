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


def sum_S2rates_over_iso(rates_dict):
    """
    rates_dict is a nested dictionary with the structure:
    {
        "mix1": {
            "iso1": [val1, val2, ...],
            "iso2": [...],
            ...
        },
        "mix2": {
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
    result = {}

    for mix, iso_dict in rates_dict.items():
        sum_vector = None
        for values in iso_dict.values():
            if sum_vector is None:
                sum_vector = [0.0] * len(values)
            for i, val in enumerate(values):
                sum_vector[i] += val
        result[mix] = sum_vector

    return result


def define_OECD_NEA_PHASE_IIIB_mixes():
    """
    Define individual mixes for OECD NEA Phase IIIB benchmark.
    Each mix is a dictionary with isotopes and their respective concentrations.
    Returns a dictionary with mix numbers as keys and their compositions as values.
    The isotopes and their concentrations are based on the OECD NEA Phase IIIB benchmark specifications.
    MIXES_def.keys() are the "ROD_ID_NUMBER" defined in the OECD NEA Phase IIIB benchmark specifications.
    """


    MIXES_def = {
                1:{"U234": 1.0443E-05, "U235": 1.1284E-03, "U236": 6.9317E-06, "U238": 2.1606E-02, "O16": 4.5504E-02},
                2:{"U234": 7.5720E-06, "U235": 8.2904E-04, "U236": 5.1701E-06, "U238": 2.1907E-02, "O16": 4.5497E-02},
                3:{"U234": 6.2468E-06, "U235": 6.9087E-04, "U236": 4.3570E-06, "U238": 2.2046E-02, "O16": 4.5494E-02},
                4:{"U234": 4.7008E-06, "U235": 5.2968E-04, "U236": 3.4083E-06, "U238": 2.2208E-02, "O16": 4.5491E-02},
                5:{"U234": 5.8824E-06, "U235": 6.5057E-04, "U236": 4.1028E-06, "U238": 2.0759E-02, "O16": 4.5095E-02, "Gd154": 3.2253E-05, "Gd155": 2.2141E-04, "Gd156": 3.0778E-04, "Gd157": 2.3576E-04, "Gd158": 3.7393E-04, "Gd160": 3.3200E-04},
                }
    return MIXES_def

