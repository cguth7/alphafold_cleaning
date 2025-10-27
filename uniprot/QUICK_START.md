# Quick Start - Run This on Your Local Machine

**The environment where this was built blocks UniProt API access. You need to run this on your local machine.**

---

## Why You Need This

Your existing PubTator data shows papers that **mention** genes.

This UniProt pipeline gets papers that are **ABOUT** genes (curated annotations).

**Big difference:**
- PubTator: "This paper mentions TP53 in context of disease X"
- UniProt: "This paper characterizes TP53 function, structure, or mechanism"

---

## Run It Locally (3 Commands)

```bash
# 1. Get the code
git pull origin claude/uniprot-citation-extraction-011CUXokXRDHhh5g4fKfKyu2
cd uniprot

# 2. Test if UniProt is accessible
python test_api_access.py
# Should show: "✓ UniProt API is accessible!"

# 3. Run the full pipeline
./run_pipeline.sh

# Done! Check data/outputs/ for your files
```

---

## What You'll Get

```
data/outputs/
├── monthly_publications_long.csv          # Monthly papers per protein
├── protein_summaries.csv                  # Stats per protein
├── yearly_publications.csv                # Yearly aggregates
└── global_monthly_timeseries.csv          # Overall trends

Expected:
- 50,000-200,000 publications
- ~20,000 proteins
- Full history (1900s-present)
- ~1-2 hour runtime
```

---

## Key Output File

**`monthly_publications_long.csv`**

```csv
accession,gene_name,protein_name,year_month,total_pubs,curated_pubs,mapped_pubs
P04637,TP53,Cellular tumor antigen p53,2020-01,45,12,33
P04637,TP53,Cellular tumor antigen p53,2020-02,38,8,30
P04637,TP53,Cellular tumor antigen p53,2020-03,52,15,37
...
```

This tells you **exactly which papers are about which proteins**, not just mentioning them.

---

## Compare to Your Existing Data

### Your PubTator Data (2020-2023)
- 19,690 proteins
- Disease-relevant papers that mention genes
- Good for: "Research involving gene X in disease context"

### New UniProt Data (Full History)
- ~20,000 proteins
- Papers that characterize/study genes
- Good for: "Research focused on gene X specifically"

---

## Questions?

Read the full documentation: `README.md`

Or just run `./run_pipeline.sh` and it'll work! 🚀
