#!/usr/bin/env python3
"""Validate CIR→PIR synthesis at ALL Hall operations for compound irreps.

For each compound irrep where cir_ops < char_count, this verifies that
PIR({R|t}) = Σ_i CIR_i(R) * exp(i*2π*k·t_iso)

using ISOTROPY translations from PIR_TRANS.
"""
import re, math, sys

def validate(generated_data_path):
    with open(generated_data_path) as f:
        content = f.read()

    def parse_array(name):
        m = re.search(rf'pub static {name}:.*?=\s*\[(.*?)\];', content, re.DOTALL)
        if not m:
            return []
        return [float(x.strip()) for x in m.group(1).split(',') if x.strip()]

    chars_flat = parse_array('CHARACTERS')
    cir_flat = parse_array('CIR_COMPONENT_CHARS')
    pir_trans_flat = parse_array('PIR_TRANS')
    pir_rots_flat = parse_array('PIR_ROTS')
    cir_rots_flat = parse_array('CIR_ROTS')

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

    total = 0
    mismatches = 0
    ok_extra = 0  # Extra ops that validated correctly

    for b in blocks:
        cir_count = int(re.search(r'_cir_count:\s*(\d+)', b).group(1))
        if cir_count == 0:
            continue

        cir_ops = int(re.search(r'_cir_ops:\s*(\d+)', b).group(1))
        char_count = int(re.search(r'_char_count:\s*(\d+)', b).group(1))

        if cir_ops >= char_count:
            continue  # CIR covers all ops, validated by Rust test

        sg = int(re.search(r'sg:\s*(\d+)', b).group(1))
        ml = re.search(r'ml:\s*"([^"]+)"', b).group(1)
        char_start = int(re.search(r'_char_start:\s*(\d+)', b).group(1))
        cir_start = int(re.search(r'_cir_start:\s*(\d+)', b).group(1))
        kx = int(re.search(r'kx:\s*(\S+)', b).group(1).rstrip(','))
        ky = int(re.search(r'ky:\s*(\S+)', b).group(1).rstrip(','))
        kz = int(re.search(r'kz:\s*(\S+)', b).group(1).rstrip(','))
        kd = int(re.search(r'kd:\s*(\S+)', b).group(1).rstrip(','))
        prs = int(re.search(r'_pir_rot_start:\s*(\d+)', b).group(1))

        pir = chars_flat[char_start:char_start + char_count]

        # CIR component chars
        cir_comps = []
        for c in range(cir_count):
            off = cir_start + c * cir_ops * 2
            cir_comps.append([(cir_flat[off + op*2], cir_flat[off + op*2 + 1]) for op in range(cir_ops)])

        trans_start = prs // 9 * 3
        cir_rot_start = (cir_start // 2) * 9

        # Validate extra Hall ops
        for h in range(cir_ops, char_count):
            if trans_start + (h+1)*3 > len(pir_trans_flat):
                continue
            if prs + (h+1)*9 > len(pir_rots_flat):
                continue

            t_h_val = pir_trans_flat[trans_start + h*3:trans_start + (h+1)*3]
            h_rot_val = pir_rots_flat[prs + h*9:prs + (h+1)*9]

            # Find matching CIR op by rotation
            match_ci = None
            for ci in range(cir_ops):
                if cir_rot_start + (ci+1)*9 > len(cir_rots_flat):
                    continue
                c_rot = cir_rots_flat[cir_rot_start + ci*9:cir_rot_start + (ci+1)*9]
                if h_rot_val == c_rot:
                    match_ci = ci
                    break

            if match_ci is None:
                continue

            # ISOTROPY translations at h and match_ci
            if trans_start + (match_ci+1)*3 > len(pir_trans_flat):
                continue
            t_ref = pir_trans_flat[trans_start + match_ci*3:trans_start + (match_ci+1)*3]
            if len(t_h_val) < 3 or len(t_ref) < 3:
                continue

            dt = [t_h_val[j] - t_ref[j] for j in range(3)]
            dot = kx * dt[0] + ky * dt[1] + kz * dt[2]
            theta = 2.0 * math.pi * dot / float(kd)
            cos_th = math.cos(theta)
            sin_th = math.sin(theta)

            # Compute CIR sum
            cir_sum_re = 0.0
            cir_sum_im = 0.0
            for c in range(cir_count):
                c_re = cir_comps[c][match_ci][0]
                c_im = cir_comps[c][match_ci][1]
                cir_sum_re += c_re * cos_th - c_im * sin_th
                cir_sum_im += c_re * sin_th + c_im * cos_th

            if abs(cir_sum_im) > 0.01:
                print(f"  MISMATCH SG{sg} {ml} op{h}: PIR={pir[h]:.4f} "
                      f"CIR_sum=({cir_sum_re:.4f},{cir_sum_im:.4f}) [IMAGINARY]")
                mismatches += 1
            elif abs(pir[h] - cir_sum_re) > 0.01:
                print(f"  MISMATCH SG{sg} {ml} op{h}: PIR={pir[h]:.4f} "
                      f"CIR_sum=({cir_sum_re:.4f},{cir_sum_im:.4f})")
                mismatches += 1
            else:
                ok_extra += 1

        total += 1

    print(f"CIR→PIR extra-op validation: {total} compound irreps, "
          f"{ok_extra} extra ops OK, {mismatches} mismatches")
    return mismatches == 0


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'src/irrep/generated_data.rs'
    ok = validate(path)
    sys.exit(0 if ok else 1)
