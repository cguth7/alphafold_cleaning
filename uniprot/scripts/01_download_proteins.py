#!/usr/bin/env python3
"""
Download all human protein entries from UniProt.

This script downloads both reviewed (Swiss-Prot) and unreviewed (TrEMBL)
human protein entries from the UniProt REST API using streaming endpoints
for efficiency.

Usage:
    python 01_download_proteins.py [--reviewed-only] [--output-dir DATA_DIR]

Author: AlphaFold Impact Analysis Pipeline
Date: October 2025
"""

import argparse
import logging
import requests
import json
import time
from pathlib import Path
from datetime import datetime
from typing import Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_proteins.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# UniProt REST API base URL
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"

# Human taxonomy ID
HUMAN_TAXON_ID = "9606"


def download_with_retry(url: str, output_file: Path, max_retries: int = 4) -> bool:
    """
    Download data from UniProt with exponential backoff retry logic.

    Args:
        url: API endpoint URL
        output_file: Path to save downloaded data
        max_retries: Maximum number of retry attempts

    Returns:
        bool: True if successful, False otherwise
    """
    retry_delays = [2, 4, 8, 16]  # Exponential backoff in seconds

    # Add proper headers as required by UniProt API
    headers = {
        'User-Agent': 'AlphaFold-Pipeline/1.0 (research@example.com)',
        'Accept': 'application/json'
    }

    for attempt in range(max_retries + 1):
        try:
            logger.info(f"Downloading from {url} (attempt {attempt + 1}/{max_retries + 1})")

            # Stream the response to handle large files
            response = requests.get(url, headers=headers, stream=True, timeout=300)
            response.raise_for_status()

            # Write to file in chunks
            with open(output_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            # Verify file was written
            if output_file.exists() and output_file.stat().st_size > 0:
                logger.info(f"Successfully downloaded to {output_file} ({output_file.stat().st_size:,} bytes)")
                return True
            else:
                logger.error(f"Downloaded file is empty: {output_file}")
                return False

        except requests.exceptions.RequestException as e:
            logger.warning(f"Download attempt {attempt + 1} failed: {e}")

            if attempt < max_retries:
                delay = retry_delays[attempt]
                logger.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                logger.error(f"All {max_retries + 1} download attempts failed")
                return False

        except Exception as e:
            logger.error(f"Unexpected error during download: {e}")
            return False

    return False


def download_human_proteins(output_dir: Path, reviewed_only: bool = False) -> dict:
    """
    Download human protein data from UniProt.

    Args:
        output_dir: Directory to save downloaded files
        reviewed_only: If True, only download reviewed (Swiss-Prot) entries

    Returns:
        dict: Summary statistics of the download
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    stats = {
        'start_time': datetime.now().isoformat(),
        'reviewed_downloaded': False,
        'unreviewed_downloaded': False,
        'reviewed_file': None,
        'unreviewed_file': None,
        'reviewed_size_bytes': 0,
        'unreviewed_size_bytes': 0
    }

    # Download reviewed proteins (Swiss-Prot)
    logger.info("="*60)
    logger.info("Downloading REVIEWED human proteins (Swiss-Prot)...")
    logger.info("="*60)

    reviewed_url = (
        f"{UNIPROT_BASE_URL}/stream?"
        f"query=organism_id:{HUMAN_TAXON_ID}%20AND%20reviewed:true&"
        f"format=json"
    )
    reviewed_file = output_dir / "human_reviewed.json"

    if download_with_retry(reviewed_url, reviewed_file):
        stats['reviewed_downloaded'] = True
        stats['reviewed_file'] = str(reviewed_file)
        stats['reviewed_size_bytes'] = reviewed_file.stat().st_size

        # Quick count of entries
        try:
            with open(reviewed_file, 'r') as f:
                data = json.load(f)
                reviewed_count = len(data.get('results', []))
                stats['reviewed_count'] = reviewed_count
                logger.info(f"Reviewed proteins: {reviewed_count:,} entries")
        except Exception as e:
            logger.warning(f"Could not count reviewed entries: {e}")

    # Download unreviewed proteins (TrEMBL) unless reviewed_only flag is set
    if not reviewed_only:
        logger.info("="*60)
        logger.info("Downloading UNREVIEWED human proteins (TrEMBL)...")
        logger.info("="*60)
        logger.info("Note: This may take longer as it's a much larger dataset")

        unreviewed_url = (
            f"{UNIPROT_BASE_URL}/stream?"
            f"query=organism_id:{HUMAN_TAXON_ID}%20AND%20reviewed:false&"
            f"format=json"
        )
        unreviewed_file = output_dir / "human_unreviewed.json"

        if download_with_retry(unreviewed_url, unreviewed_file):
            stats['unreviewed_downloaded'] = True
            stats['unreviewed_file'] = str(unreviewed_file)
            stats['unreviewed_size_bytes'] = unreviewed_file.stat().st_size

            # Quick count of entries
            try:
                with open(unreviewed_file, 'r') as f:
                    data = json.load(f)
                    unreviewed_count = len(data.get('results', []))
                    stats['unreviewed_count'] = unreviewed_count
                    logger.info(f"Unreviewed proteins: {unreviewed_count:,} entries")
            except Exception as e:
                logger.warning(f"Could not count unreviewed entries: {e}")

    stats['end_time'] = datetime.now().isoformat()

    # Save download stats
    stats_file = output_dir / "download_stats.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Download statistics saved to {stats_file}")

    return stats


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Download human protein entries from UniProt"
    )
    parser.add_argument(
        '--reviewed-only',
        action='store_true',
        help='Only download reviewed (Swiss-Prot) entries'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/raw'),
        help='Output directory for downloaded files (default: data/raw)'
    )

    args = parser.parse_args()

    logger.info("="*60)
    logger.info("UniProt Human Proteins Download")
    logger.info("="*60)
    logger.info(f"Output directory: {args.output_dir.absolute()}")
    logger.info(f"Reviewed only: {args.reviewed_only}")
    logger.info("")

    # Execute download
    stats = download_human_proteins(args.output_dir, args.reviewed_only)

    # Print summary
    logger.info("="*60)
    logger.info("DOWNLOAD SUMMARY")
    logger.info("="*60)
    logger.info(f"Reviewed proteins: {'✓' if stats['reviewed_downloaded'] else '✗'}")
    if stats['reviewed_downloaded']:
        logger.info(f"  File: {stats['reviewed_file']}")
        logger.info(f"  Size: {stats['reviewed_size_bytes']:,} bytes ({stats['reviewed_size_bytes']/1024/1024:.1f} MB)")
        if 'reviewed_count' in stats:
            logger.info(f"  Count: {stats['reviewed_count']:,} entries")

    if not args.reviewed_only:
        logger.info(f"Unreviewed proteins: {'✓' if stats['unreviewed_downloaded'] else '✗'}")
        if stats['unreviewed_downloaded']:
            logger.info(f"  File: {stats['unreviewed_file']}")
            logger.info(f"  Size: {stats['unreviewed_size_bytes']:,} bytes ({stats['unreviewed_size_bytes']/1024/1024:.1f} MB)")
            if 'unreviewed_count' in stats:
                logger.info(f"  Count: {stats['unreviewed_count']:,} entries")

    logger.info("="*60)
    logger.info("Download complete!")
    logger.info("="*60)


if __name__ == "__main__":
    main()
