#!/usr/bin/env python3
"""Check ISOTROPY PIR translations vs spglib Hall translations."""
import json, re, math

with open("src/irrep/generated_data.rs") as f:
    content = f.read()

with open("scripts/hall_operations.json") as f:
    hall_data = json.load(f)

# Build sg_num → list of (hall_num, rots, trans)
sg_halls = {}
for hall_str, hd in hall_data.items():
    sg = hd["sg"]
    if sg not in sg_halls:
        sg_halls[sg] = []
    sg_halls[sg].append((int(hall_str), hd["rots"], hd["trans"]))

def parse_array(name):
    m = re.search(rf'pub static {name}:.*?=\s*\[(.*?)\];', content, re.DOTALL)
    return [float(x.strip()) for x in m.group(1).split(',') if x.strip()]

pir_trans = parse_array('PIR_TRANS')
pir_rots = parse_array('PIR_ROTS')

blocks = []
depth = 0; current = ''
for line in content.split('\n'):
    if 'IrrepRecord {' in line: current = line + '\n'; depth = 1
    elif depth > 0:
        current += line + '\n'
        depth += line.count('{') - line.count('}')
        if depth == 0: blocks.append(current)

# For each SG's first irrep (GM1 typically), compare PIR trans with Hall trans
sg_done = set()
total_match = 0
total_mismatch = 0
mismatch_examples = []

for b in blocks:
    sg = int(re.search(r'sg:\s*(\d+)', b).group(1))
    if sg in sg_done:
        continue
    sg_done.add(sg)

    prs = int(re.search(r'_pir_rot_start:\s*(\d+)', b).group(1))
    cc = int(re.search(r'_char_count:\s*(\d+)', b).group(1))
    ml = re.search(r'ml:\s*"([^"]+)"', b).group(1)

    # Get PIR rotations and translations
    trans_start = prs // 9 * 3
    if trans_start + cc * 3 > len(pir_trans):
        continue

    # Find matching Hall setting
    if sg not in sg_halls:
        continue

    best_match = None
    best_count = 0
    for hall_num, hall_rots, hall_trans in sg_halls[sg]:
        if len(hall_rots) != cc:
            continue
        # Check if rotations match
        match_count = 0
        for op in range(cc):
            p_rot = [int(x) for x in pir_rots[prs+op*9:prs+(op+1)*9]]
            h_rot = hall_rots[op]
            if p_rot == h_rot:
                match_count += 1
        if match_count == cc:
            best_match = (hall_num, hall_trans)
            best_count = match_count
            break

    if best_match is None:
        continue

    hall_num, hall_trans = best_match

    # Compare translations
    diffs = []
    for op in range(cc):
        p_t = pir_trans[trans_start+op*3:trans_start+(op+1)*3]
        h_t = hall_trans[op]
        dt = [abs(p_t[j] - h_t[j]) for j in range(3)]
        max_dt = max(dt)
        if max_dt > 0.001:
            diffs.append((op, p_t, h_t, max_dt))

    if diffs:
        total_mismatch += 1
        if len(mismatch_examples) < 10:
            mismatch_examples.append((sg, ml, cc, diffs[:3]))
    else:
        total_match += 1

print(f"ISOTROPY vs spglib translation check:")
print(f"  {total_match} SGs: translations match")
print(f"  {total_mismatch} SGs: translations differ")
print(f"  {230 - total_match - total_mismatch} SGs: no matching Hall setting found")

if mismatch_examples:
    print(f"\nMismatch examples:")
    for sg, ml, cc, diffs in mismatch_examples:
        print(f"  SG{sg} {ml} ({cc} ops):")
        for op, p_t, h_t, dt in diffs[:3]:
            print(f"    op{op}: PIR=[{p_t[0]:.4f},{p_t[1]:.4f},{p_t[2]:.4f}] "
                  f"Hall=[{h_t[0]:.4f},{h_t[1]:.4f},{h_t[2]:.4f}] diff={dt:.4f}")
