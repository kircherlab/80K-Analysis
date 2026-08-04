### SCREEN v4 overlap and enrichment
- download the screen v4 bed files from the SCREEN portal and filter them with the command below
```bash
cd ./data/SCREEN_tissues

OUTDIR="enrichment_over_tissues_filtered"
mkdir -p "$OUTDIR"

for f in *.noccl.cCREs.bed.gz; do
  tissue="${f%.noccl.cCREs.bed.gz}"
  out="${OUTDIR}/${tissue}_cCRE_open_filtered.bed.gz"

  zcat "$f" \
	| awk '$0 ~ /(CA-H3K4me3|CA-only|CA-TF|dELS|pELS|PLS)/' \
    | gzip -c > "$out"

  echo "Wrote: $out";
done
```
- the raw files from the screen v4 were too big but the filtered files are provided.