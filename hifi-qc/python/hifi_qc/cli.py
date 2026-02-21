import click


@click.command()
@click.argument("bam_files", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("-o", "--output", default="hifi_qc_report.html", help="Output HTML report path")
@click.option("-r", "--reference", default=None, type=click.Path(exists=True), help="Reference FASTA for GC bias")
@click.option("-t", "--threads", default=None, type=int, help="Number of threads")
@click.option("--sample-fraction", default=None, type=float, help="Fraction of reads to sample (0.0-1.0)")
@click.option("--seed", default=None, type=int, help="RNG seed for reproducible sampling")
@click.option("--sample-names", default=None, type=str, help="Comma-separated sample names")
def main(bam_files, output, reference, threads, sample_fraction, seed, sample_names):
    """Quality control for PacBio HiFi sequencing data."""
    from hifi_qc._core import process_bam_files

    results = process_bam_files(
        list(bam_files),
        reference=reference,
        threads=threads,
        sample_fraction=sample_fraction,
        seed=seed,
    )
    click.echo(f"Processed {len(bam_files)} file(s). Results: {len(results)} sample(s).")
    click.echo(f"Report would be written to: {output}")
