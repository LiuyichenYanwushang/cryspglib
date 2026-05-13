#!/usr/bin/env python3
"""Debug the reorder mapping for SG9 M1M2."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Monkey-patch generate_irrep_data to print debug for SG9 M1M2
import generate_irrep_data as gen

# Re-run the data generation with debug
orig_reorder = gen._reorder_to_spglib_order

def debug_reorder(sg, ml, chars_flat, char_starts, char_counts,
                  matrices_flat, mat_starts, mat_counts,
                  pir_rots_flat, pir_rot_starts, rots_map,
                  spinor_irreps, spinor_starts, spinor_counts,
                  cir_comp_flat=None, cir_comp_rots=None,
                  cir_comp_starts=None, cir_comp_counts=None, cir_comp_ops=None,
                  kvec_map=None):
    """Wrapper that prints debug for SG9 M1M2."""
    for i in range(len(ml)):
        if sg[i] == 9 and ml[i] == 'M1M2':
            n_ops = char_counts[i]
            pir_rots = rots_map.get((9, 'M1M2'), [])
            print(f"DEBUG SG9 M1M2: n_ops={n_ops} len(pir_rots)={len(pir_rots)}")
            if pir_rots:
                print(f"  First ISOTROPY rot: {pir_rots[0]}")
                print(f"  All ISOTROPY rots:")
                for j, r in enumerate(pir_rots[:8]):
                    print(f"    ISO[{j}]: {r}")
            break

    result = orig_reorder(sg, ml, chars_flat, char_starts, char_counts,
                          matrices_flat, mat_starts, mat_counts,
                          pir_rots_flat, pir_rot_starts, rots_map,
                          spinor_irreps, spinor_starts, spinor_counts,
                          cir_comp_flat, cir_comp_rots,
                          cir_comp_starts, cir_comp_counts, cir_comp_ops,
                          kvec_map=kvec_map)

    # After reorder, check M1M2 PIR values
    for i in range(len(ml)):
        if sg[i] == 9 and ml[i] == 'M1M2':
            cs = char_starts[i]
            cc = char_counts[i]
            print(f"DEBUG after reorder: char_start={cs} char_count={cc}")
            for op in range(min(cc, 4)):
                print(f"  PIR[{op}] = {chars_flat[cs + op]:.4f}")
            # Check mapping
            reorder_map, hall_choice = result
            if i < len(reorder_map) and reorder_map[i] is not None:
                print(f"  mapping = {reorder_map[i]}")
            break

    return result

gen._reorder_to_spglib_order = debug_reorder

# Run generation
gen.main()
