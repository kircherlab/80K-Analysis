"""
This script removes the adapter from the design file and replaces it with a random adapter sequence: returns two different design files.

:Author: Kilian Salomon
:Contact: kilian.salomon@bih-charite.de
:Date: *29.07.2024
"""
# Example: python /home/kisa/coding/80K_MPRA/80K-Analysis/07_quality_control/scripts/remove_adapter_and_random_adapter.py --input /home/kisa/coding/80K_MPRA/design_data/design_no_duplicates_sequence_and_header.fa --random-adapter True --remove-adapter True --start 16 --length 270 --random-seed 42

import pyfastx
import click
import random
import os



# options
@click.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Fasta design file.",
)
@click.option(
    "--start",
    "start",
    required=True,
    type=int,
    help="Start of the sequence to look at (1 based; i.e. start after adapter).",
)
@click.option(
    "--length",
    "length",
    required=True,
    type=int,
    help="length of the sequence to lok at (i.e. without adapter).",
)
@click.option(
    "--random-adapter",
    "random_adapter",
    type=bool,
    default=True,
    help="Add random adapter sequences.",
)
@click.option(
    "--remove-adapter",
    "remove_adapter",
    type=bool,
    default=True,
    help="Remove adapter sequences.",
)
@click.option(
    "--random-seed",
    "random_seed",
    type=int,
    default=42,
    help="Random seed.",
)
@click.option(
    "--output",
    "output_dir",
    type=click.Path(exists=False, writable=True),
    default=".",
    help="Output directory.",
)
@click.option(
    "--filename",
    "filename",
    type=str,
    default="output",
    help="Output filename. (Default: same as input filename)",
)


def cli(input_file, start, length, random_adapter, remove_adapter, random_seed, output_dir, filename):

    def generate_random_dna(length):
        """
        Generates a random DNA sequence of a given length.

        Parameters:
            length (int): The length of the DNA sequence to generate.

        Returns:
            str: A string representing the random DNA sequence.
        """
        nucleotides = ['A', 'T', 'C', 'G']
        return ''.join(random.choice(nucleotides) for _ in range(length))

    # set random seed
    random.seed(42)

    header = []
    base_sequences = []

    # read fasta file
    fa = pyfastx.Fasta(input_file)
    # remove adapter (0-start of sequence - 2 (start is 1 based)) and (end - length of sequence)
    for i in range(len(fa)):
        header.append(fa[i].name)
        base_sequences.append(fa[i].seq[start-1:start+length-1])

    if remove_adapter:
        # write fasta file again
        output_filename = filename if filename != "output" else input_file.split(".")[0]
        output_file = output_filename.split(".")[0] + "_no_adapter.fa" if output_dir == "." else os.path.join(output_dir, output_filename.split(".")[0] + "_no_adapter.fa")
        click.echo(f"Removing adapter sequences and writing fasta to {output_file}...")
        with open(output_file, "w") as f:
            for i in range(len(header)):
                f.write(f">{header[i]}\n{base_sequences[i]}\n")

    # write fasta file with random adapter
    if not random_adapter:
        click.echo("No random adapter sequences added.")
        return

    output_filename = filename if filename != "output" else input_file.split(".")[0]
    output_file = output_filename.split(".")[0] + "_random_adapter.fa" if output_dir == "." else os.path.join(output_dir, output_filename.split(".")[0] + "_random_adapter.fa")
    click.echo(f"Adding random adapter sequences and writing fasta to {output_file}...")
    with open(output_file, "w") as f:
        for i in range(len(header)):
            adapter_start = generate_random_dna(start-1)
            adapter_end = generate_random_dna(start-1)
            f.write(f">{header[i]}\n{adapter_start}{base_sequences[i]}{adapter_end}\n")


if __name__ == "__main__":
    cli()