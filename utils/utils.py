import pandas as pd


def verbose_print(print_string, threshold, verbose_level, end='\n'):
    # print given string, if verbose_level is higher than threshold
    #
    # input print_string: str, threshold: int, verbose_level: int, end: str
    # no return

    if verbose_level > threshold:
        print(print_string, end=end)


def round_self(value, decimals):
    # simple decimal rounding function (python by itself has a tendency to round fragmented with the built-in function)
    #
    # input value: float, decimals: int
    # return rounded_value: float/int

    # If decimal less than 1, the resulting value will be an integer
    if pd.isna(value):
        return float("Nan")
    if decimals < 1:
        rounded_value = int(f"{int((value * (10 ** decimals)) + .5) / (10 ** decimals):.{decimals}f}")
        return rounded_value
    # Else, the resulting value will be a float
    else:
        rounded_value = float(f"{int((value * (10 ** decimals)) + .5) / (10 ** decimals):.{decimals}f}")
        return rounded_value
