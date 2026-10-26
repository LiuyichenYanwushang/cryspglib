#!/usr/bin/env python3
"""Check all magnetic monoclinic raw lattices/origins against the bundled iso.

Run: python3 scripts/verify_magnetic_monoclinic_frames.py
The oracle's SHOW BASIS at DISPLAY PARENT independently fixes the raw primitive
order. This is a frame check, not an embedding search or an OP-stabilizer proof.
Manual: https://iso.byu.edu/isotropy_doc.php (SETTING INTERNATIONAL, SHOW BASIS).
"""
import re
from fractions import Fraction

from generate_irrep_data import read_file, get_sections, parse_ints, parse_labels
from verify_isotropy_oracle import (
    ISO_DIR, PRIMITIVE_BASIS, run_iso_session, same_lattice,
    to_conventional, to_conventional_basis, primitive_of_conventional_cell,
)


def vectors(text):
    return [[Fraction(v) for v in row.split(',')]
            for row in text.strip('()').split('),(')]


def verify():
    raw = read_file('data_irreps.txt')
    sections = get_sections(raw)
    labels = parse_labels(raw, sections, 'irrep_label')
    parents = parse_ints(raw, sections, 'irrep_space_group')
    raw = read_file('data_magnetic.txt')
    sections = get_sections(raw)
    values = {key: parse_ints(raw, sections, key) for key in (
        'mag_iso_irrep', 'mag_iso_subgroup', 'mag_iso_basis', 'mag_iso_origin')}
    directions = parse_labels(raw, sections, 'mag_iso_orderparam_label')
    numbers = parse_labels(raw, sections, 'mag_nlabel')
    rows = {}
    for k, irrep in enumerate(values['mag_iso_irrep']):
        sg = parents[irrep - 1]
        if not 3 <= sg <= 15:
            continue
        b = values['mag_iso_basis'][9*k:9*k+9]
        o = values['mag_iso_origin'][4*k:4*k+4]
        key = (sg, labels[irrep - 1], directions[k])
        assert key not in rows, key
        rows[key] = (numbers[values['mag_iso_subgroup'][k] - 1],
                     [b[i:i+3] for i in range(0, 9, 3)],
                     [Fraction(x, o[3]) for x in o[:3]])
    seen = set()
    nonzero = 0
    for sg in range(3, 16):
        canonical = PRIMITIVE_BASIS['C' if sg in (5, 8, 9, 12, 15) else 'P']
        primitive = [canonical[2], canonical[0], canonical[1]]
        parent = run_iso_session(ISO_DIR + '/iso', [
            'PAGE 1000', 'SC 400', f'VALUE PARENT {sg}',
            'SHOW BASIS', 'DISPLAY PARENT', 'QUIT'])
        parent_rows = [vectors(line.strip().strip('*'))
                       for line in parent.splitlines()
                       if line.strip().strip('*').startswith('(')]
        assert parent_rows == [primitive], (sg, parent_rows, primitive)
        text = run_iso_session(ISO_DIR + '/iso', [
            'PAGE 1000', 'SC 400', 'SET I ALL OR 1', 'SETTING MAGNETIC',
            f'VALUE PARENT {sg}', 'SHOW IRREP', 'SHOW SUBGROUP', 'SHOW BASIS',
            'SHOW ORIGIN', 'SHOW SIZE', 'SHOW DIRECTION', 'DISPLAY ISOTROPY', 'QUIT'])
        for line in text.splitlines():
            match = re.fullmatch(
                r'(\S+)\s+(\d+\.\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)',
                line.strip().strip('*'))
            if match is None:
                continue
            ml, number, symbol, _, direction, basis, origin = match.groups()
            key = (sg, ml, direction)
            assert key in rows and key not in seen, key
            seen.add(key)
            expected_number, b, o = rows[key]
            assert number == expected_number, key
            expected = to_conventional_basis(b, primitive)
            actual = primitive_of_conventional_cell(vectors(basis), PRIMITIVE_BASIS[symbol[0]])
            assert same_lattice(expected, actual), (key, expected, actual)
            expected_origin = to_conventional(o, primitive)
            assert expected_origin == vectors(origin)[0], (key, expected_origin, origin)
            nonzero += any(expected_origin)
    assert seen == rows.keys(), rows.keys() - seen
    assert (len(seen), nonzero) == (262, 105)
    print('262 official monoclinic rows; 105 nonzero origins; 13 parent frames; 0 mismatches')


if __name__ == '__main__':
    verify()
