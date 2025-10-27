#!/usr/bin/env python3
"""
Quick test of UniProt API access and data download for a small sample.

This script tests if UniProt API is accessible and downloads data for
just 10 well-known proteins as a proof of concept.

If this works, the full pipeline will work.
"""

import requests
import json
import sys

# Test proteins (well-known genes)
TEST_PROTEINS = [
    'P04637',  # TP53
    'P01308',  # INS
    'P05231',  # IL6
    'P01375',  # TNF
    'P01137',  # TGFB1
    'P00533',  # EGFR
    'P04626',  # ERBB2
    'P35354',  # PTGS2
    'P42574',  # CASP3
    'P31749',  # AKT1
]

def test_single_protein(accession):
    """Test downloading a single protein."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    headers = {
        'User-Agent': 'AlphaFold-Pipeline/1.0',
        'Accept': 'application/json'
    }

    try:
        response = requests.get(url, headers=headers, timeout=10)
        if response.status_code == 200:
            data = response.json()
            # Extract basic info
            gene_name = data.get('genes', [{}])[0].get('geneName', {}).get('value', 'Unknown')
            num_refs = len(data.get('references', []))
            return {
                'accession': accession,
                'gene_name': gene_name,
                'num_references': num_refs,
                'status': 'SUCCESS'
            }
        else:
            return {
                'accession': accession,
                'status': f'FAILED: {response.status_code}'
            }
    except Exception as e:
        return {
            'accession': accession,
            'status': f'ERROR: {str(e)}'
        }

def main():
    print("="*60)
    print("UniProt API Access Test")
    print("="*60)
    print(f"Testing {len(TEST_PROTEINS)} well-known proteins...\n")

    results = []
    for acc in TEST_PROTEINS:
        print(f"Testing {acc}...", end=" ")
        result = test_single_protein(acc)
        results.append(result)

        if result['status'] == 'SUCCESS':
            print(f"✓ {result['gene_name']}: {result['num_references']} references")
        else:
            print(f"✗ {result['status']}")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    successful = sum(1 for r in results if r['status'] == 'SUCCESS')
    print(f"Successful: {successful}/{len(TEST_PROTEINS)}")

    if successful == len(TEST_PROTEINS):
        print("\n✓ UniProt API is accessible!")
        print("✓ Full pipeline should work")

        # Save sample data
        with open('test_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nSample data saved to test_results.json")
        return 0
    else:
        print("\n✗ UniProt API is NOT accessible")
        print("✗ Need to run on different network")
        return 1

if __name__ == '__main__':
    sys.exit(main())
