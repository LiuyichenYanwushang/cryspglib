#!/bin/bash
# Full irrep data regeneration pipeline.
#
# Steps:
#   1. Rust:  export hall operations → hall_operations.json (if outdated)
#   2. Python: parse ISO/PIR/CIR zips + reorder to spglib Hall order → generated_data.rs
#
# Usage: ./scripts/regenerate_all.sh
# After running: cargo test --package cryspglib

set -euo pipefail
cd "$(dirname "$0")/.."

echo "=== Step 1: Export hall operations ==="
cargo test --package cryspglib --lib irrep::corep::tests::export_hall_operations -- --nocapture 2>&1 | tail -1
echo "  done"

echo "=== Step 2: Generate irrep data (with spglib reorder + CIR padding) ==="
python3 scripts/generate_irrep_data.py
echo "  done"

echo ""
echo "=== Pipeline complete ==="
echo "Run: cargo test --package cryspglib"
