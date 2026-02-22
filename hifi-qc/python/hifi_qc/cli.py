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
    from hifi_qc.report import generate_report

    names = sample_names.split(",") if sample_names else None

    click.echo(f"Processing {len(bam_files)} file(s)...")
    results = process_bam_files(
        list(bam_files),
        reference=reference,
        threads=threads,
        sample_fraction=sample_fraction,
        seed=seed,
        sample_names=names,
    )

    click.echo(f"Generating report for {len(results)} sample(s)...")
    generate_report(results, output)
    click.echo(f"Report written to: {output}")
