"""Mapping from isotropy direction codes to order-parameter component strings.

The string is only a *label* for the direction: the authoritative selector is
the ISOTROPY direction label (`P1`, `C2`, `S1`, …), which this repository always
carries as `IsotropyRecord::direction_label`.

Provenance of every entry below:

* `dim = 1` and `dim = 3` component patterns were read from the bundled program with
  `SHOW SUBGROUP` + `SHOW DIRECTION VECTOR` + `DISPLAY ISOTROPY` (command in
  `scripts/verify_isotropy_oracle.py`), so they are the program's own strings.
  They were re-checked over the `dim = 3` records of a 175-irrep sample spanning
  all crystal systems (`target/scratch/verify_dim3_map.py`): 781 label checks,
  0 mismatches for `P1/P2/P3/C1/S1`, and `C2` matched the cubic form
  `(a,a,b)` for cubic parents while trigonal/hexagonal parents printed the
  complex form `(a;b;a)` (22 mismatches before the split was made explicit).
  The program uses `;` where the components are complex.  The Rust selector
  and oracle compare these patterns ignoring whitespace and treating `;` as
  `,`; the stored strings do not preserve the program's separators everywhere.
* `dim = 2` entries are **cryspglib's internal** component notation; the
  official components vary per irrep (e.g. SG 194 GM6+ P2 is `(a,0.577a)`).
* `dim >= 4` entries are **cryspglib's own** compact notation
  (`LABEL(free)/DIMD`), not a string the archive stores: the program prints full
  component lists (e.g. `(a,0,b,0,-b,0)`) that are per-irrep, so a table keyed
  by `(dim, free, label)` cannot reproduce them.  Users who need the program's
  exact spelling should select by `Label`.
"""

# (dim, free, label) -> component string, for the cases where the program's
# string is a function of the triple alone.  Verified as described above.
OFFICIAL = {
    (1, 1, "P1"): "(a)",
    (3, 1, "P1"): "(a,0,0)",
    (3, 1, "P2"): "(a,a,0)",
    (3, 1, "P3"): "(a,a,a)",
    (3, 2, "C1"): "(a,b,0)",
    (3, 2, "C2"): "(a,a,b)",  # cubic parents; see NON_CUBIC_C2
    (3, 3, "S1"): "(a,b,c)",
}

# Trigonal/hexagonal parents print the complex form instead.
NON_CUBIC_C2 = "(a;b;a)"


def _crystal_system_is_cubic(sg_number):
    return sg_number >= 195


def direction_str(dim, free, label, parent_sg=None):
    """Component string for one isotropy record (see the module docstring)."""
    official = OFFICIAL.get((dim, free, label))
    if official is not None:
        if (dim, free, label) == (3, 2, "C2") and parent_sg is not None:
            if not _crystal_system_is_cubic(parent_sg):
                return NON_CUBIC_C2
        return official

    if dim <= 3:
        # dim = 2 is deliberately the internal notation: the program's strings
        # depend on whether the parent irrep is complex (SG 5 `L1` prints
        # `(a;a)` while SG 91 `A1` prints `(a,0)`), so no triple-keyed table is
        # correct.  These strings stay stable for `Descriptor` matching.
        if dim == 1:
            return "(a)"
        if dim == 2:
            if free == 1:
                return {"P1": "(a,0)", "P2": "(0,a)", "P3": "(a,a)", "P4": "(a,-a)"}.get(
                    label, f"(P:{label})"
                )
            return "(a,b)"
        return f"(dim{dim})"

    # dim >= 4: cryspglib's compact notation, never claimed to be the program's.
    return f"{label}({free})/{dim}D"


def build_direction_map(dir_vals, dim_vals, free_vals, label_vals, parent_sgs=None):
    """Map each direction code to a description string.

    Convenience wrapper over [`direction_str`] for callers that group by code.
    Prefer `direction_str` per record: a code can be shared by cubic and
    non-cubic parents for the family-dependent `C2` spelling.
    """
    code_map = {}
    for i in range(len(dir_vals)):
        code = dir_vals[i]
        if code in code_map:
            continue
        parent_sg = parent_sgs[i] if parent_sgs is not None and i < len(parent_sgs) else None
        code_map[code] = direction_str(
            dim_vals[i] if i < len(dim_vals) else 0,
            free_vals[i] if i < len(free_vals) else 0,
            label_vals[i].strip() if i < len(label_vals) else "?",
            parent_sg,
        )
    return code_map
