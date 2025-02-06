# How to?
> in order to postprocess the mismatches table of an Allo-Count run.

Execute `multi_launcher.sh` which will run one of the 3 .py scripts in this directory based on the filter type wanted.

Example: for **bed** filtering:
*Note: requires the Conda environment activated.*

```bash
./multi_launcher.sh bed ../tutorial/example.bed output ../output/runs/run/run_tables/
```

Arguments are:
1. filter type: "bed", "rsID", "genes-transcripts"
2. filter file: examples files are [there](https://github.com/huguesrichard/Allopipe/tree/main/tutorial)
3. output directory
4. mismatches table path directory