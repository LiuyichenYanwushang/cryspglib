#!/usr/bin/env python3
"""Check Bloch phase between ISOTROPY and spglib for non-Gamma k-points."""
import json, re, math

with open("src/irrep/generated_data.rs") as f:
    content = f.read()
with open("scripts/hall_operations.json") as f:
    hall_data = json.load(f)

sg_halls = {}
for hall_str, hd in hall_data.items():
    sg = hd["sg"]
    if sg not in sg_halls: sg_halls[sg] = []
    sg_halls[sg].append((int(hall_str), hd["rots"], hd["trans"]))

def parse_array(name):
    m = re.search(rf'pub static {name}:.*?=\s*\[(.*?)\];', content, re.DOTALL)
    return [float(x.strip()) for x in m.group(1).split(',') if x.strip()]

pir_trans = parse_array('PIR_TRANS')
pir_rots = parse_array('PIR_ROTS')
chars = parse_array('CHARACTERS')

blocks = []; depth = 0; current = ''
for line in content.split('\n'):
    if 'IrrepRecord {' in line: current = line + '\n'; depth = 1
    elif depth > 0:
        current += line + '\n'
        depth += line.count('{') - line.count('}')
        if depth == 0: blocks.append(current)

# Count phases
real_phases = 0
imag_phases = 0
examples = []

for b in blocks:
    sg = int(re.search(r'sg:\s*(\d+)', b).group(1))
    kx = int(re.search(r'kx:\s*(\S+)', b).group(1).rstrip(','))
    ky = int(re.search(r'ky:\s*(\S+)', b).group(1).rstrip(','))
    kz = int(re.search(r'kz:\s*(\S+)', b).group(1).rstrip(','))
    kd = int(re.search(r'kd:\s*(\S+)', b).group(1).rstrip(','))

    if kx == 0 and ky == 0 and kz == 0:
        continue  # Skip Gamma point

    prs_m = re.search(r'_pir_rot_start:\s*(\d+)', b)
    cc_m = re.search(r'_char_count:\s*(\d+)', b)
    cs_m = re.search(r'_char_start:\s*(\d+)', b)
    ml_m = re.search(r'ml:\s*"([^"]+)"', b)
    if not prs_m or not cc_m or not cs_m or not ml_m:
        continue

    prs = int(prs_m.group(1))
    cc = int(cc_m.group(1))
    cs = int(cs_m.group(1))
    ml = ml_m.group(1)

    if sg not in sg_halls:
        continue

    # Find matching Hall setting
    hall_trans = None
    for hall_num, hall_rots, h_trans in sg_halls[sg]:
        if len(hall_rots) != cc:
            continue
        ok = True
        for op_i in range(cc):
            p_rot = [int(x) for x in pir_rots[prs+op_i*9:prs+(op_i+1)*9]]
            if p_rot != hall_rots[op_i]:
                ok = False; break
        if ok:
            hall_trans = h_trans; break
    if hall_trans is None:
        continue

    trans_start = prs // 9 * 3
    for op_i in range(cc):
        p_t = pir_trans[trans_start+op_i*3:trans_start+(op_i+1)*3]
        h_t = hall_trans[op_i]
        dt = [h_t[j] - p_t[j] for j in range(3)]
        if max(abs(d) for d in dt) < 0.001:
            continue
        dot = kx*dt[0] + ky*dt[1] + kz*dt[2]
        theta = 2*math.pi*dot/kd
        phase_cos = math.cos(theta)
        phase_sin = math.sin(theta)
        if abs(phase_sin) > 0.01:
            imag_phases += 1
            if len(examples) < 5:
                examples.append((sg, ml, op_i, kx, ky, kz, kd, p_t, h_t, dt, theta, chars[cs+op_i]))
        else:
            real_phases += 1

print(f"Non-Gamma operations with translation diff:")
print(f"  {real_phases} have real phase (cos=+-1)")
print(f"  {imag_phases} have imaginary phase (problem!)")
print()

if examples:
    print("Imaginary phase examples:")
    for sg, ml, op_i, kx, ky, kz, kd, p_t, h_t, dt, theta, pir_val in examples:
        print(f"  SG{sg} {ml} op{op_i}: k=({kx},{ky},{kz})/{kd}")
        print(f"    PIR_t={[f'{t:.4f}' for t in p_t]} Hall_t={[f'{t:.4f}' for t in h_t]}")
        print(f"    dt={[f'{d:.4f}' for d in dt]} theta={theta:.4f} phase_cos={math.cos(theta):.4f} phase_sin={math.sin(theta):.4f}")
        print(f"    PIR_char={pir_val:.4f}")
