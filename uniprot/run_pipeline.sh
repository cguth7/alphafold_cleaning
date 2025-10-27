#!/bin/bash
# Master script to run the entire UniProt citation extraction pipeline
# Usage: ./run_pipeline.sh [--reviewed-only]

set -e  # Exit on error

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
REVIEWED_ONLY=""
if [ "$1" == "--reviewed-only" ]; then
    REVIEWED_ONLY="--reviewed-only"
    echo -e "${YELLOW}Running in REVIEWED ONLY mode (Swiss-Prot only)${NC}"
fi

# Function to print section headers
print_header() {
    echo ""
    echo -e "${BLUE}================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================================================${NC}"
    echo ""
}

# Function to print step completion
print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

# Function to print errors
print_error() {
    echo -e "${RED}✗ $1${NC}"
}

# Start pipeline
print_header "UniProt Citation Extraction Pipeline"
echo "Start time: $(date)"
START_TIME=$(date +%s)

# Check Python
if ! command -v python &> /dev/null; then
    print_error "Python not found. Please install Python 3.8 or higher."
    exit 1
fi

PYTHON_VERSION=$(python --version 2>&1 | awk '{print $2}')
echo "Python version: $PYTHON_VERSION"

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo -e "${YELLOW}Virtual environment not found. Creating one...${NC}"
    python -m venv venv
    print_success "Virtual environment created"
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate || source venv/Scripts/activate

# Install requirements
echo "Checking dependencies..."
pip install -q -r requirements.txt
print_success "Dependencies installed"

# Step 1: Download
print_header "STEP 1/4: Downloading UniProt Data"
python scripts/01_download_proteins.py $REVIEWED_ONLY
if [ $? -eq 0 ]; then
    print_success "Download completed"
else
    print_error "Download failed"
    exit 1
fi

# Step 2: Parse
print_header "STEP 2/4: Parsing Publications"
python scripts/02_parse_publications.py
if [ $? -eq 0 ]; then
    print_success "Parsing completed"
else
    print_error "Parsing failed"
    exit 1
fi

# Step 3: Aggregate
print_header "STEP 3/4: Aggregating Monthly Data"
python scripts/03_aggregate_monthly.py
if [ $? -eq 0 ]; then
    print_success "Aggregation completed"
else
    print_error "Aggregation failed"
    exit 1
fi

# Step 4: Validate
print_header "STEP 4/4: Validating Data Quality"
python scripts/04_validate_data.py
if [ $? -eq 0 ]; then
    print_success "Validation passed"
else
    print_error "Validation failed (but pipeline completed)"
fi

# Calculate runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
MINUTES=$((RUNTIME / 60))
SECONDS=$((RUNTIME % 60))

# Final summary
print_header "PIPELINE COMPLETE"
echo "End time: $(date)"
echo "Total runtime: ${MINUTES}m ${SECONDS}s"
echo ""
echo "Output files are in:"
echo "  - data/raw/           (downloaded JSON files)"
echo "  - data/processed/     (parsed publications)"
echo "  - data/outputs/       (aggregated datasets)"
echo ""
echo "Next steps:"
echo "  1. Review validation results: cat data/outputs/validation_results.json"
echo "  2. Explore data: see notebooks/exploratory_analysis.ipynb"
echo "  3. Read documentation: README.md"
echo ""
print_success "All done!"
