import click

from .talon_to_counts import write_top_tables
from .talon_validate import (
    datasets_present,
    min_le_max,
    novelty_sets_ok,
    qc_identity_ok,
    schema_ok,
    spotcheck_reads,
)


@click.group()
def cli():
    "se4s-isoform CLI"


@cli.command()
@click.option("--tsv", required=True)
@click.option("--qc", required=True)
@click.option("--datasets", default="ENCFF003OWX,ENCFF019HRC,ENCFF669LWV,ENCFF676BYQ")
@click.option(
    "--enable-spotcheck",
    is_flag=True,
    default=False,
    help="Enable SAM round-trip spot check (disabled by default). Requires --se4s-root",
)
@click.option("--se4s-root", default="")
def validate(tsv, qc, datasets, enable_spotcheck, se4s_root):
    """Run a quick validation suite on a TALON read-wise TSV and QC log.

    By default the SAM round-trip spotcheck is disabled (it requires access to
    the `work/longread/bulk/<dataset>` directories). Use `--enable-spotcheck`
    and provide `--se4s-root` to enable that check.
    """
    expected = set(datasets.split(","))
    checks = {
        "schema_ok": schema_ok(tsv),
        "coords_min<=max": min_le_max(tsv),
        "datasets_present": datasets_present(tsv, expected),
        "novelty_sets_ok": novelty_sets_ok(tsv),
        "qc_identity_ok": qc_identity_ok(qc),
    }
    if enable_spotcheck:
        if not se4s_root:
            raise click.BadParameter(
                "--se4s-root must be set when --enable-spotcheck is used"
            )
        checks["spotcheck_reads"] = spotcheck_reads(tsv, se4s_root)

    for k, v in checks.items():
        click.echo(f"{k}: {v}")
    click.echo("VALID" if all(checks.values()) else "INVALID")


@cli.command()
@click.option("--tsv", required=True)
@click.option("--out", required=True, help="Output directory for CSVs")
@click.option("--dataset", default=None)
@click.option("--top", type=int, default=50)
def counts(tsv, out, dataset, top):
    tx_csv, gene_csv = write_top_tables(tsv, out, dataset=dataset, top=top)
    click.echo(f"Wrote {tx_csv}")
    click.echo(f"Wrote {gene_csv}")
