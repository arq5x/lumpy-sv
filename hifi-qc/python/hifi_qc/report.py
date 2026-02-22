"""HTML report generator for hifi-qc."""

import json
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

COLORS = [
    "#e41a1c",  # red
    "#377eb8",  # blue
    "#4daf4a",  # green
    "#984ea3",  # purple
    "#ff7f00",  # orange
    "#a65628",  # brown
    "#f781bf",  # pink
    "#999999",  # gray
]


def generate_report(samples: list[dict], output_path: str) -> None:
    """Generate a self-contained HTML QC report from sample metrics.

    Parameters
    ----------
    samples : list[dict]
        List of sample metric dicts as returned by ``process_bam_files``.
    output_path : str
        Destination file path for the rendered HTML report.
    """
    template_dir = Path(__file__).parent / "templates"
    env = Environment(loader=FileSystemLoader(str(template_dir)), autoescape=False)
    template = env.get_template("report.html.j2")

    for i, sample in enumerate(samples):
        sample["color"] = COLORS[i % len(COLORS)]

    html = template.render(
        samples=samples,
        samples_json=json.dumps(samples, default=_json_serializer),
        colors=COLORS,
        has_coverage=any("genome_mean_coverage" in s for s in samples),
        has_gc_bias=any("gc_bias_coverage_by_gc" in s for s in samples),
    )
    Path(output_path).write_text(html)


def _json_serializer(obj):
    """Fallback serializer for numpy-like numeric types."""
    if hasattr(obj, "__float__"):
        return float(obj)
    if hasattr(obj, "__int__"):
        return int(obj)
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
