"""This is our main entry point"""
import logging

import click
from Bio.SeqIO import read

from .nnk import ReferenceLibrary, analyze_library
from .util import get_verbosity_level


@click.group("nnk")
@click.pass_context
def nnk(ctx):
    pass


@nnk.command("reference")
@click.pass_context
@click.option(
    "-v", "--verbose", count=True, default=5, help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option("--primer-begin", required=True, help="The primer begin postion")
@click.option("--primer-end", required=True, help="The primer end postion")
@click.option(
    "--start-codon",
    required=True,
    help="Relative to the end of the start of the primer, where does the reading frame start",
    type=int,
)
@click.option(
    "--end-codon",
    required=True,
    help="Relative to the start of the end primer, where does the reading frame end",
    type=int,
)
@click.option(
    "--reference",
    required=True,
    help="The reference sequence to use that the library was created. Should be in a fasta fle. Will be stripped down to reading frame from primers and codons",
    type=click.Path(file_okay=True, dir_okay=False, exists=True),
)
@click.argument(
    "input", required=True, type=click.Path(file_okay=True, dir_okay=False, writable=True, resolve_path=True)
)
@click.argument(
    "output",
    required=False,
    type=click.Path(file_okay=True, dir_okay=False, writable=True, resolve_path=True),
    default="library.json.gz",
)
def reference(ctx, verbose, primer_begin, primer_end, start_codon, end_codon, reference, input, output):
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    click.echo(
        "Creating reference library with primer {}{}-->{}{}".format(
            primer_begin, start_codon * "N", end_codon * -1 * "N", primer_end
        )
    )
    reference_seq = str(read(reference, "fasta").seq)
    reference_library = ReferenceLibrary(primer_begin, primer_end, reference_seq, start_codon, end_codon)
    reference_library.apply(input)
    reference_library.to_json(output)


@nnk.command("analyze")
@click.pass_context
@click.option(
    "-v", "--verbose", count=True, default=5, help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
    "--reference-lib",
    required=True,
    help="The reference library to use from the reference library generator",
    type=click.Path(file_okay=True, dir_okay=False, exists=True),
)
@click.argument("input", required=True, type=click.Path(file_okay=True, dir_okay=False, exists=True))
@click.argument("output", required=True, type=click.Path(file_okay=True, dir_okay=False, exists=False))
def analyze(ctx, verbose, reference_lib, input, output):
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)
    reference = ReferenceLibrary.from_json(reference_lib)
    library = analyze_library(reference, input)
    library.to_csv(output)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    analyze()
    # reference()
