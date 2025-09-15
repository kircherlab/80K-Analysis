# score_variants.py

from io import StringIO
from alphagenome import colab_utils
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
import pandas as pd
from tqdm import tqdm
import os # Import os for checking file existence

ALPHAGENOME_API_KEY = "your_api_key_here"  # Replace with your actual API key
dna_model = dna_client.create(ALPHAGENOME_API_KEY)

# Load VCF file containing variants.
# This path needs to be accessible from where the Slurm job runs.
vcf_file = '/sc-projects/sc-proj-bih-reg-seqs/users/kisa11/projects/alphaGenome/all_80k_variants_vcf_100.tsv'
vcf_file = '/sc-projects/sc-proj-bih-reg-seqs/users/kisa11/projects/alphaGenome/all_80k_variants_vcf.tsv'
output_dir = '/sc-projects/sc-proj-bih-reg-seqs/users/kisa11/projects/alphaGenome/results'

# Check if the VCF file exists
if not os.path.exists(vcf_file):
    raise FileNotFoundError(f"VCF file not found at: {vcf_file}. Please ensure the path is correct and accessible.")

vcf = pd.read_csv(vcf_file, sep='\t')

required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
for column in required_columns:
    if column not in vcf.columns:
        raise ValueError(f'VCF file is missing required column: {column}.')

organism = 'human'  # @param ["human", "mouse"] {type:"string"}

# @markdown Specify length of sequence around variants to predict:
sequence_length = '1MB'  # @param ["2KB", "16KB", "100KB", "500KB", "1MB"] { type:"string" }
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
    f'SEQUENCE_LENGTH_{sequence_length}'
]

# Scorer selections
score_rna_seq = False
score_cage = False
score_procap = False
score_atac = True
score_dnase = False
score_chip_histone = False
score_chip_tf = False
score_polyadenylation = False
score_splice_sites = False
score_splice_site_usage = False
score_splice_junctions = False

# Parse organism specification.
organism_map = {
    'human': dna_client.Organism.HOMO_SAPIENS,
    'mouse': dna_client.Organism.MUS_MUSCULUS,
}
organism = organism_map[organism]

# Parse scorer specification.
scorer_selections = {
    'rna_seq': score_rna_seq,
    'cage': score_cage,
    'procap': score_procap,
    'atac': score_atac,
    'dnase': score_dnase,
    'chip_histone': score_chip_histone,
    'chip_tf': score_chip_tf,
    'polyadenylation': score_polyadenylation,
    'splice_sites': score_splice_sites,
    'splice_site_usage': score_splice_site_usage,
    'splice_junctions': score_splice_junctions,
}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]

# Remove any scorers or output types that are not supported for the chosen organism.
unsupported_scorers = [
    scorer
    for scorer in selected_scorers
    if (
        organism.value
        not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
    )
    | (
        (scorer.requested_output == dna_client.OutputType.PROCAP)
        & (organism == dna_client.Organism.MUS_MUSCULUS)
    )
]
if len(unsupported_scorers) > 0:
  print(
      f'Excluding {unsupported_scorers} scorers as they are not supported for'
      f' {organism}.'
  )
  for unsupported_scorer in unsupported_scorers:
    selected_scorers.remove(unsupported_scorer)


# Score variants in the VCF file.
chunk_size = 100
results = []
chunk_num = 0

for i, vcf_row in tqdm(vcf.iterrows(), total=len(vcf)):

    variant = genome.Variant(
        chromosome=str(vcf_row.CHROM),
        position=int(vcf_row.POS),
        reference_bases=vcf_row.REF,
        alternate_bases=vcf_row.ALT,
        name=vcf_row.variant_id,
    )

    interval = variant.reference_interval.resize(sequence_length)

    variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=selected_scorers,
        organism=organism,
    )

    results.append(variant_scores)

    if (i + 1) % chunk_size == 0:
        df_chunk = variant_scorers.tidy_scores(results)

        chunk_filename = f"{output_dir}/chunk_{chunk_num:05d}_{len(results)}_variants_80k_alphaGenome_score.tsv.gz"
        df_chunk.to_csv(chunk_filename, sep="\t", index=False)
        print(f"Saved {len(results)} variants to {chunk_filename}")

        results = []  # Clear memory
        chunk_num += 1

# Write remaining variants
if results:
    df_chunk = variant_scorers.tidy_scores(results)
    chunk_filename = f"{output_dir}/chunk_{chunk_num:05d}_{len(results)}_variants_80k_alphaGenome_score.tsv.gz"
    df_chunk.to_csv(chunk_filename, sep="\t", index=False)
    print(f"Saved {len(results)} variants to {chunk_filename}")





















