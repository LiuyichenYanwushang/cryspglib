#!/usr/bin/env python3
"""Debug SG9 M1M2 CIR→PIR phase."""
import re, math, json, zipfile

with open("src/irrep/generated_data.rs") as f:
    content = f.read()

def parse_array(name):
    m = re.search(rf'pub static {name}:.*?=\s*\[(.*?)\];', content, re.DOTALL)
    return [float(x.strip()) for x in m.group(1).split(',') if x.strip()]

chars = parse_array('CHARACTERS')
cir = parse_array('CIR_COMPONENT_CHARS')
pir_trans = parse_array('PIR_TRANS')
pir_rots = parse_array('PIR_ROTS')
cir_rots = parse_array('CIR_ROTS')

# Find SG9 M1M2
blocks = []; depth = 0; current = ''
for line in content.split('\n'):
    if 'IrrepRecord {' in line: current = line + '\n'; depth = 1
    elif depth > 0:
        current += line + '\n'
        depth += line.count('{') - line.count('}')
        if depth == 0: blocks.append(current)

for b in blocks:
    if 'sg: 9,' not in b or 'M1M2' not in b:
        continue
    cc = int(re.search(r'_char_count:\s*(\d+)', b).group(1))
    cs = int(re.search(r'_char_start:\s*(\d+)', b).group(1))
    co = int(re.search(r'_cir_ops:\s*(\d+)', b).group(1))
    cstart = int(re.search(r'_cir_start:\s*(\d+)', b).group(1))
    ccount = int(re.search(r'_cir_count:\s*(\d+)', b).group(1))
    prs = int(re.search(r'_pir_rot_start:\s*(\d+)', b).group(1))
    kx = int(re.search(r'kx:\s*(\S+)', b).group(1).rstrip(','))
    ky = int(re.search(r'ky:\s*(\S+)', b).group(1).rstrip(','))
    kz = int(re.search(r'kz:\s*(\S+)', b).group(1).rstrip(','))
    kd = int(re.search(r'kd:\s*(\S+)', b).group(1).rstrip(','))

    print(f"SG9 M1M2: char_count={cc} cir_ops={co} cir_count={ccount}")
    print(f"  k=({kx},{ky},{kz})/{kd}")

    # PIR chars
    print(f"  PIR chars: {[chars[cs+op] for op in range(cc)]}")

    # PIR translations (ISOTROPY)
    trans_start = prs // 9 * 3
    print(f"  PIR ISOTROPY trans (from index {trans_start}):")
    for op in range(cc):
        t = pir_trans[trans_start+op*3:trans_start+(op+1)*3]
        print(f"    op{op}: [{t[0]:.4f},{t[1]:.4f},{t[2]:.4f}]")

    # PIR rotations
    print(f"  PIR ISOTROPY rots:")
    for op in range(cc):
        r = pir_rots[prs+op*9:prs+(op+1)*9]
        print(f"    op{op}: {r}")

    # CIR chars and rotations
    cir_rot_start = (cstart // 2) * 9
    print(f"  CIR (first {co} ops):")
    for op in range(co):
        print(f"    op{op}: chars=", end="")
        for comp in range(ccount):
            off = cstart + comp*co*2 + op*2
            print(f"comp{comp}=({cir[off]:.2f},{cir[off+1]:.2f}) ", end="")
        r = cir_rots[cir_rot_start+op*9:cir_rot_start+(op+1)*9]
        t = pir_trans[trans_start+op*3:trans_start+(op+1)*3]
        print(f"rot={r} t=[{t[0]:.4f},{t[1]:.4f},{t[2]:.4f}]")

    # For each extra Hall op, find the matching CIR op
    print(f"\n  Validation for extra ops ({co}..{cc-1}):")
    for h in range(co, cc):
        h_rot = pir_rots[prs + h*9:prs + (h+1)*9]
        t_h = pir_trans[trans_start + h*3:trans_start + (h+1)*3]

        # Find matching CIR
        match_ci = None
        for ci in range(co):
            c_rot = cir_rots[cir_rot_start + ci*9:cir_rot_start + (ci+1)*9]
            if h_rot == c_rot:
                match_ci = ci
                break

        print(f"    op{h}: rot={h_rot} t_iso=[{t_h[0]:.4f},{t_h[1]:.4f},{t_h[2]:.4f}]")
        print(f"      match_ci={match_ci}")

        if match_ci is not None:
            # Get t_ref from the matching CIR position
            t_ref = pir_trans[trans_start + match_ci*3:trans_start + (match_ci+1)*3]
            dt = [t_h[j] - t_ref[j] for j in range(3)]
            dot = kx*dt[0] + ky*dt[1] + kz*dt[2]
            theta = 2*math.pi*dot/float(kd)
            print(f"      t_ref=[{t_ref[0]:.4f},{t_ref[1]:.4f},{t_ref[2]:.4f}]")
            print(f"      dt={dt} k·dt={dot:.4f} θ={theta:.4f} cos={math.cos(theta):.4f} sin={math.sin(theta):.4f}")

            # Compute CIR sum at this op
            cir_sum_re = 0.0; cir_sum_im = 0.0
            for comp in range(ccount):
                off = cstart + comp*co*2 + match_ci*2
                c_re = cir[off]; c_im = cir[off+1]
                new_re = c_re*math.cos(theta) - c_im*math.sin(theta)
                new_im = c_re*math.sin(theta) + c_im*math.cos(theta)
                print(f"      CIR comp{comp}: ({c_re:.2f},{c_im:.2f}) * phase = ({new_re:.4f},{new_im:.4f})")
                cir_sum_re += new_re; cir_sum_im += new_im
            print(f"      CIR_sum=({cir_sum_re:.4f},{cir_sum_im:.4f}) vs PIR={chars[cs+h]:.4f}")

            # Also compare PIR at match_ci and h
            print(f"      PIR[op{match_ci}]={chars[cs+match_ci]:.4f} PIR[op{h}]={chars[cs+h]:.4f} ratio={chars[cs+h]/chars[cs+match_ci]:.4f}")

    break

# Also check: what are the original ISOTROPY PIR data for M1M2?
with zipfile.ZipFile("isotropy_subgroup/PIR_data.zip") as zf:
    with zf.open("PIR_data.txt") as f:
        text = f.read().decode('utf-8', errors='replace')

for i, line in enumerate(text.split('\n')):
    if '9' in line.split()[:3] and 'M1M2' in line:
        parts = line.split()
        print(f"\n  ISOTROPY PIR source M1M2: k=({parts[4]},{parts[5]},{parts[6]})/{parts[7]}")
        rots_line = list(map(int, text.split('\n')[i+1].split()))
        trans_line = list(map(float, text.split('\n')[i+2].split()))
        chars_line = list(map(float, text.split('\n')[i+3].split()))
        n = len(rots_line)//9
        print(f"  n_ops={n}")
        for op in range(min(6, n)):
            r = rots_line[op*9:(op+1)*9]
            t = trans_line[op*3:(op+1)*3]
            print(f"    ISO op{op}: rot={r[:6]}... t=[{t[0]:.4f},{t[1]:.4f},{t[2]:.4f}] char={chars_line[op]:.4f}")
        break
