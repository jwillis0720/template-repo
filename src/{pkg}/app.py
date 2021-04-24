"""This is our main entry point"""
import logging
import sys
import click

@click.group()
@click.pass_context
def {pkg{(ctx):
    pass


@sadie.command()
@click.pass_context
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=5,
    help="Vebosity level, ex. -vvvvv for debug level logging",
)
@click.option(
)
@click.argument(
    "input",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True, readable=True, resolve_path=True),
)
@click.argument(
    "output", required=False, type=click.Path(file_okay=True, dir_okay=False, writable=True, resolve_path=True)
)
def cmd1(ctx, verbose, input,output):
    numeric_level = get_verbosity_level(verbose)
    logging.basicConfig(level=numeric_level)

if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    cmd1()
