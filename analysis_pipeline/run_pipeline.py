#!/usr/bin/env python3
"""
run_pipeline.py - Master Pipeline for AlphaFold2 Impact Analysis

This script runs the complete analysis pipeline:
1. CEM Matching - Creates matched treatment/control groups
2. Semester Aggregation - Aggregates to semester level, creates ΔY variables
3. Event Study - Runs difference-in-differences event studies
4. DOI Integration - Prepares DOI data for Shi & Evans metrics

Usage:
    python run_pipeline.py           # Run full pipeline
    python run_pipeline.py --step 1  # Run only step 1 (CEM matching)
    python run_pipeline.py --step 3  # Run only step 3 (event studies)

Author: Replicating Danilo Messinese's Stata analysis in Python
"""

import sys
import argparse
from pathlib import Path
import time

# Add scripts directory to path
SCRIPT_DIR = Path(__file__).parent / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))


def run_step(step_num, step_name, module_name):
    """Run a pipeline step with timing."""
    print("\n" + "=" * 70)
    print(f"STEP {step_num}: {step_name}")
    print("=" * 70 + "\n")

    start = time.time()

    try:
        module = __import__(module_name)
        module.main()
        elapsed = time.time() - start
        print(f"\n✓ Step {step_num} completed in {elapsed:.1f} seconds")
        return True
    except Exception as e:
        print(f"\n✗ Step {step_num} failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description="AlphaFold2 Impact Analysis Pipeline")
    parser.add_argument('--step', type=int, help='Run only a specific step (1-4)')
    parser.add_argument('--skip-doi', action='store_true', help='Skip DOI integration step')
    args = parser.parse_args()

    steps = [
        (1, "CEM Matching", "01_cem_matching"),
        (2, "Semester Aggregation", "02_semester_aggregation"),
        (3, "Event Study Analysis", "03_event_study"),
        (4, "DOI Integration", "04_integrate_dois"),
    ]

    print("\n" + "#" * 70)
    print("#" + " " * 68 + "#")
    print("#" + "    ALPHAFOLD2 IMPACT ANALYSIS PIPELINE".center(68) + "#")
    print("#" + " " * 68 + "#")
    print("#" + "    Replicating Danilo's CEM + Event Study Analysis".center(68) + "#")
    print("#" + " " * 68 + "#")
    print("#" * 70)

    total_start = time.time()
    results = []

    for step_num, step_name, module_name in steps:
        # Skip if specific step requested
        if args.step and args.step != step_num:
            continue

        # Skip DOI step if requested
        if args.skip_doi and step_num == 4:
            continue

        success = run_step(step_num, step_name, module_name)
        results.append((step_num, step_name, success))

        if not success and not args.step:
            print(f"\nPipeline stopped at step {step_num}")
            break

    # Summary
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)

    for step_num, step_name, success in results:
        status = "✓ SUCCESS" if success else "✗ FAILED"
        print(f"  Step {step_num}: {step_name:30s} {status}")

    total_elapsed = time.time() - total_start
    print(f"\nTotal time: {total_elapsed:.1f} seconds")

    # Output locations
    print("\n" + "-" * 70)
    print("OUTPUT LOCATIONS")
    print("-" * 70)
    base = Path(__file__).parent
    print(f"  Data:    {base / 'data'}")
    print(f"  Outputs: {base / 'outputs'}")
    print(f"  Figures: {base / 'figures'}")

    print("\n" + "#" * 70)
    print("# PIPELINE COMPLETE".ljust(69) + "#")
    print("#" * 70 + "\n")


if __name__ == "__main__":
    main()
