#!/bin/bash
# Full irrep data regeneration pipeline.
#
# Steps:
#   1. Python: parse ISO/PIR/CIR zips → generated_data.rs (without reorder data)
#   2. Rust:  generate char_reorder_map.txt (ISOTROPY→spglib mapping)
#   3. Python: apply_char_reorder.py → reorder_data.rs
#   4. Append reorder_data.rs to generated_data.rs
#
# Usage: ./scripts/regenerate_all.sh
# After running: cargo test --package cryspglib --lib irrep

set -euo pipefail
cd "$(dirname "$0")/.."

echo "=== Step 1: Generate main data ==="
python3 scripts/generate_irrep_data.py
echo "  done"

echo "=== Step 2: Generate char reorder map ==="
cargo test --package cryspglib --lib irrep::corep::tests::generate_char_reorder_map -- --nocapture 2>&1 | tail -1
echo "  done"

echo "=== Step 3: Build reorder data ==="
python3 scripts/apply_char_reorder.py
echo "  done"

echo "=== Step 4: Append reorder data to generated_data.rs ==="
# Save original (without reorder) and restore
python3 -c "
import subprocess
result = subprocess.run(['git', 'show', 'HEAD:src/irrep/generated_data.rs'],
                       capture_output=True, text=True)
with open('src/irrep/generated_data.rs', 'w') as f:
    f.write(result.stdout)
with open('src/irrep/reorder_data.rs') as f:
    reorder = f.read()
with open('src/irrep/generated_data.rs', 'a') as f:
    f.write('\n')
    f.write(reorder)
print('  appended reorder data')
"
echo "  done"

echo ""
echo "=== Pipeline complete ==="
echo "Run: cargo test --package cryspglib --lib irrep"
