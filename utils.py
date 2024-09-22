def smart_round(coeff, tolerance=1e-10):
    """Round to integer or keep 9 decimal."""
    if abs(coeff - round(coeff)) < tolerance:
        return round(coeff)
    else:
        return round(coeff, 10)


def get_sign(coeff, is_first_term):
    if is_first_term:
        return "-" if coeff < 0 else ""
    return "-" if coeff < 0 else "+"


def get_format_coeff(coeff, power):
    abs_coeff = abs(coeff)
    if abs_coeff == 1 and power != 0:
        coeff_str = ""
    else:
        coeff_str = f"{abs_coeff}"
    return coeff_str


def get_format_power(coeff_str, power):
    if power == 0:
        return f"{coeff_str}"
    elif power == 1:
        return f"{coeff_str} * x" if coeff_str else "x"
    else:
        return f"{coeff_str} * x**{power}" if coeff_str else f"x**{power}"


def ensemble_polynomial(an, use_smart_round=True, keep_all_zeros=False, keep_first_zeros=True):
    """
    Ensemble px by input coefficient an, make it clean and beauty.

    """
    terms = []
    n = len(an) - 1
    skip_zero = False if keep_first_zeros else True

    if use_smart_round:
        an = [smart_round(a) for a in an]

    for i in range(n, -1, -1):
        coeff = an[i]
        is_first_term = i == n

        if keep_all_zeros:
            pass
        elif skip_zero == bool(coeff):
            pass
        elif coeff == 0:
            continue
        else:
            skip_zero = True

        sign = get_sign(coeff, is_first_term)
        coeff_str = get_format_coeff(coeff, i)
        term = get_format_power(coeff_str, i)
        terms.append(f" {sign} {term}" if sign else f"{term}")
    px = "".join(terms)
    return px


if __name__ == "__main__":
    an = [1.9999999999999425, -2.4999999999999996, 1]
    px = ensemble_polynomial(an)
    print(px)

    an = [1, 0.0, 2.5, -1, 0, 1, 0, 0]
    px = ensemble_polynomial(an)
    print(px)
