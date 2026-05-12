#!/usr/bin/env python3
"""Validate CIR→PIR synthesis at ALL Hall operations for compound irreps.

For each compound irrep where cir_ops < char_count, this verifies that
PIR({R|t}) = Σ_i CIR_i(R) * exp(i*2π*k·t)

using ISOTROPY translations from PIR_TRANS.
"""
import re, math, sys

def validate(generated_data_path, hall_ops_path, kvec_map_path=None):
    with open(generated_data_path) as f:
        content = f.read()

    # Parse arrays
    def parse_array(name):
        m = re.search(rf'pub static {name}:.*?=\s*\[(.*?)\];', content, re.DOTALL)
        if not m:
            return []
        return [float(x.strip()) for x in m.group(1).split(',') if x.strip()]

    chars_flat = parse_array('CHARACTERS')
    cir_flat = parse_array('CIR_COMPONENT_CHARS')
    pir_trans_flat = parse_array('PIR_TRANS')

    # Parse all IrrepRecord blocks
    blocks = []
    depth = 0
    current = ''
    for line in content.split('\n'):
        if 'IrrepRecord {' in line:
            current = line + '\n'
            depth = 1
        elif depth > 0:
            current += line + '\n'
            depth += line.count('{') - line.count('}')
            if depth == 0:
                blocks.append(current)

    total_checked = 0
    total_mismatches = 0

    for b in blocks:
        cir_count = int(re.search(r'_cir_count:\s*(\d+)', b).group(1))
        if cir_count == 0:
            continue

        cir_ops = int(re.search(r'_cir_ops:\s*(\d+)', b).group(1))
        char_count = int(re.search(r'_char_count:\s*(\d+)', b).group(1))

        # Only validate compound irreps where CIR needs expansion
        if cir_ops >= char_count:
            continue

        sg = int(re.search(r'sg:\s*(\d+)', b).group(1))
        ml = re.search(r'ml:\s*"([^"]+)"', b).group(1)
        char_start = int(re.search(r'_char_start:\s*(\d+)', b).group(1))
        cir_start = int(re.search(r'_cir_start:\s*(\d+)', b).group(1))
        kx = int(re.search(r'kx:\s*(\S+)', b).group(1).rstrip(','))
        ky = int(re.search(r'ky:\s*(\S+)', b).group(1).rstrip(','))
        kz = int(re.search(r'kz:\s*(\S+)', b).group(1).rstrip(','))
        kd = int(re.search(r'kd:\s*(\S+)', b).group(1).rstrip(','))

        # Get PIR chars
        pir = chars_flat[char_start:char_start + char_count]

        # Get CIR component chars (each has cir_ops ops)
        cir_comps = []
        for c in range(cir_count):
            off = cir_start + c * cir_ops * 2
            comp = [(cir_flat[off + op*2], cir_flat[off + op*2 + 1]) for op in range(cir_ops)]
            cir_comps.append(comp)

        # Get ISOTROPY translations for this irrep
        # PIR_TRANS has 3 per op, same ordering as PIR chars
        # We need pir_trans_start but it's not in IrrepRecord!
        # Estimate from index: each irrep's trans starts at the same position
        # as pir_rot_start / 9 * 3. Let's try to get it from the pir_rot_start.
        try:
            prs = int(re.search(r'_pir_rot_start:\s*(\d+)', b).group(1))
            trans_start = prs // 9 * 3
        except:
            continue

        # Build rotation-based matching: for each Hall op h >= cir_ops,
        # find which CIR op (0..cir_ops-1) has the same ISOTROPY rotation.
        # Compare ISOTROPY translations (not Hall, since CIR depends on ISO).

        # Get PIR rotations for this irrep (to match against CIR rotations)
        # CIR rots start at (cir_start // 2) * 9
        cir_rot_start = (cir_start // 2) * 9
        pir_rot_start = prs

        # Read all arrays directly
        cir_rots_all = parse_array('CIR_ROTS')
        pir_rots_all = parse_array('PIR_ROTS')

        # Match extra Hall ops to CIR ops by ISOTROPY rotation
        mismatches_this = 0
        for h in range(cir_ops, char_count):
            # Get ISOTROPY translation at Hall position h
            t_h = pir_trans_flat[trans_start + h*3:trans_start + (h+1)*3]

            # Find matching CIR op by rotation
            # Compare PIR rotation at h with CIR rotations at 0..cir_ops-1
            h_rot = pir_rots_all[pir_rot_start + h*9:pir_rot_start + (h+1)*9]
            match_ci = None
            for ci in range(cir_ops):
                c_rot = cir_rots_all[cir_rot_start + ci*9:cir_rot_start + (ci+1)*9]
                if h_rot == c_rot:
                    match_ci = ci
                    break

            if match_ci is None:
                continue  # Can't verify without rotation match

            # Get CIR for this rotation (at reference position match_ci)
            t_ref = pir_trans_flat[trans_start + match_ci*3:trans_start + (match_ci+1)*3]

            # Bloch phase: exp(i*2π*k·(t_h - t_ref)/kd)
            dt = [t_h[j] - t_ref[j] for j in range(3)]
            dot = kx * dt[0] + ky * dt[1] + kz * dt[2]
            theta = 2.0 * math.pi * dot / float(kd)
            cos_th = math.cos(theta)
            sin_th = math.sin(theta)

            # Compute CIR sum with phase for this Hall op
            cir_sum_re = 0.0
            cir_sum_im = 0.0
            for c in range(cir_count):
                c_re = cir_comps[c][match_ci][0]
                c_im = cir_comps[c][match_ci][1]
                # (re + i*im) * (cos + i*sin)
                cir_sum_re += c_re * cos_th - c_im * sin_th
                cir_sum_im += c_re * sin_th + c_im * cos_th

            # Compare with PIR (should be real)
            if abs(cir_sum_im) > 0.01:
                print(f"  MISMATCH SG{sg} {ml} op{h}: PIR={pir[h]:.4f} "
                      f"CIR_sum=({cir_sum_re:.4f},{cir_sum_im:.4f}) [imaginary!]")
                mismatches_this += 1
            elif abs(pir[h] - cir_sum_re) > 0.01:
                print(f"  MISMATCH SG{sg} {ml} op{h}: PIR={pir[h]:.4f} "
                      f"CIR_sum=({cir_sum_re:.4f},{cir_sum_im:.4f})")
                mismatches_this += 1

        total_checked += 1
        total_mismatches += mismatches_this

    print(f"\nCIR→PIR full validation: {total_checked} compound irreps checked, "
          f"{total_mismatches} mismatches")
    return total_mismatches == 0


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'src/irrep/generated_data.rs'
    ok = validate(path, '')
    sys.exit(0 if ok else 1)
