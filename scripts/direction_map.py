"""
Mapping from isotropy direction codes to human-readable direction strings.

The ISOTROPY direction codes are sequential indices that encode specific
directions in the order-parameter space. The mapping is determined by:

- Order parameter dimension (D): the dimensionality of the irrep space
- Free parameter count (f): dimensionality of the isotropy subspace
- Label (P1, C2, S1, 4D1, etc.): the direction type within D-dimensional space

Data source: isotropy_subgroup/iso_data/data_isotropy.txt
Reference: Stokes & Hatch (1988), Isotropy Subgroups of the 230
           Crystallographic Space Groups

For dim <= 3 we provide explicit component descriptions.
For dim >= 4 we use a compact format: <label>(<free>)/<dim>D
"""


def build_direction_map(dir_vals, dim_vals, free_vals, label_vals):
    """Build a dict mapping each direction code to a descriptive string.

    The input arrays are parallel (position N corresponds across all arrays).
    dir_vals contains the numeric direction codes that appear in the data.

    Args:
        dir_vals: list of direction code numbers (from isotropy_direction)
        dim_vals: list of order parameter dimensions
        free_vals: list of free parameter counts
        label_vals: list of direction labels ("P1", "C2", etc.)

    Returns:
        dict[int -> str]: direction code -> description string
    """
    # Build a unique mapping: for each direction code value, record (dim, free, label)
    # Use the first occurrence of each code to determine the mapping
    code_map = {}
    for i in range(len(dir_vals)):
        code = dir_vals[i]
        if code not in code_map:
            dim = dim_vals[i] if i < len(dim_vals) else 0
            free = free_vals[i] if i < len(free_vals) else 0
            label = label_vals[i].strip() if i < len(label_vals) else "?"
            code_map[code] = (dim, free, label)

    result = {}
    for code in sorted(code_map.keys()):
        dim, free, label = code_map[code]
        result[code] = _direction_str(dim, free, label)

    return result


def _direction_str(dim, free, label):
    """Convert (dim, free, label) to a direction description string."""
    # For dimensions 1-3, provide explicit component notation
    if dim <= 3:
        return _explicit_direction(dim, free, label)
    else:
        return _compact_direction(dim, free, label)


def _explicit_direction(dim, free, label):
    """Explicit component notation for dim <= 3."""
    if dim == 1:
        if label == "P1":
            return "(a)"
        return "(a)"  # only one possible 1D direction

    if dim == 2:
        if free == 1:
            if label == "P1":
                return "(a,0)"
            elif label == "P2":
                return "(0,a)"
            elif label == "P3":
                return "(a,a)"
            elif label == "P4":
                return "(a,-a)"
            else:
                return f"(Pn:{label})"  # unexpected P label
        elif free == 2:
            if label == "C1":
                return "(a,b)"
            else:
                return f"(a,b)"  # generic 2D
        else:
            return "(?,?)"

    if dim == 3:
        if free == 1:
            if label == "P1":
                return "(a,0,0)"
            elif label == "P2":
                return "(0,a,0)"
            elif label == "P3":
                return "(a,a,0)"
            else:
                return f"(Pn:{label})"
        elif free == 2:
            if label == "C1":
                return "(a,b,0)"
            elif label == "C2":
                return "(a,0,b)"
            else:
                return f"(Cn:{label})"
        elif free == 3:
            if label == "S1":
                return "(a,b,c)"
            else:
                return "(a,b,c)"
        else:
            return f"(dim{dim})"

    return f"(dim{dim})"


def _compact_direction(dim, free, label):
    """Compact notation for higher-dimensional directions: <label>(<free>)/<dim>D."""
    return f"{label}({free})/{dim}D"
