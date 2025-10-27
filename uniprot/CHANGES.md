# Recent Changes

## 2020-2023 Default (Latest)

**Pipeline now defaults to 2020-2023 data to match AlphaFold era analysis.**

### What Changed

**Parsing Script (`02_parse_publications.py`):**
- Added `--min-year` flag (default: 2020)
- Added `--max-year` flag (default: 2023)
- Filters publications during extraction for efficiency
- Much faster processing (~30-60 min vs ~1-2 hours)

### Why This Matters

1. **Matches your existing dataset:** PubTator data is 2020-2023
2. **Faster processing:** Less data = quicker results
3. **Direct comparison:** Same timeframe, different focus (mentions vs annotations)
4. **Still flexible:** Can override with `--min-year 1900` for full history

### Usage

**Default (2020-2023):**
```bash
./run_pipeline.sh  # Automatically uses 2020-2023
```

**Custom range:**
```bash
python scripts/02_parse_publications.py --min-year 2015 --max-year 2023
```

**Full history:**
```bash
python scripts/02_parse_publications.py --min-year 1900 --max-year 2023
```

### Expected Output

**2020-2023 (default):**
- 20,000-50,000 publications
- 48 months of data
- ~30-60 minute runtime
- Perfect for AlphaFold impact analysis

**Full history (if you want it):**
- 50,000-200,000+ publications
- 100+ years of data
- ~1-2 hour runtime
- Good for long-term baseline

### Benefits

✅ Faster extraction
✅ Matches your AlphaFold analysis period
✅ Direct comparison with PubTator data
✅ Still get gene-specific annotations (not just mentions)
✅ Can always expand date range later if needed
