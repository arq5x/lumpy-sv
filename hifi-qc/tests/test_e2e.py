"""End-to-end tests for hifi-qc."""
import os
import subprocess
import sys

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TEST_DIR, "data")
BAM_FILE = os.path.join(DATA_DIR, "test.bam")
REF_FILE = os.path.join(DATA_DIR, "test_ref.fa")
OUTPUT_DIR = os.path.join(DATA_DIR, "output")


def run_hifi_qc(*args):
    """Run hifi-qc and return the result."""
    cmd = ["hifi-qc"] + list(args)
    return subprocess.run(cmd, capture_output=True, text=True)


def test_basic_run():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output = os.path.join(OUTPUT_DIR, "basic_report.html")
    result = run_hifi_qc(BAM_FILE, "-o", output)
    assert result.returncode == 0, f"Failed: {result.stderr}"
    assert os.path.exists(output)
    with open(output) as f:
        html = f.read()
    assert "hifi-qc Report" in html
    assert "General Statistics" in html
    assert "Read Length" in html
    assert "Q-Score" in html
    assert "HiFi Passes" in html
    assert "Mapping Quality" in html
    # Coverage should be present (not sampling)
    assert "Genome" in html or "Coverage" in html
    # Per-chrom distribution section
    assert "Per-Chromosome Coverage Distribution" in html
    assert "chrom-select" in html
    assert "plot-chrom-hist" in html
    assert "plot-chrom-cumulative" in html
    # Sample filter UI
    assert "sample-search" in html
    assert "btn-select-all" in html
    assert "btn-select-none" in html
    assert "sample-checkboxes" in html
    print("PASS: basic_run")


def test_with_reference():
    output = os.path.join(OUTPUT_DIR, "ref_report.html")
    result = run_hifi_qc(BAM_FILE, "-o", output, "-r", REF_FILE)
    assert result.returncode == 0, f"Failed: {result.stderr}"
    with open(output) as f:
        html = f.read()
    assert "GC Bias" in html
    print("PASS: with_reference")


def test_sampling_mode():
    output = os.path.join(OUTPUT_DIR, "sampled_report.html")
    result = run_hifi_qc(BAM_FILE, "-o", output, "--sample-fraction", "0.5", "--seed", "42")
    assert result.returncode == 0, f"Failed: {result.stderr}"
    with open(output) as f:
        html = f.read()
    # Should still have basic sections
    assert "Read Length" in html
    print("PASS: sampling_mode")


def test_custom_sample_name():
    output = os.path.join(OUTPUT_DIR, "named_report.html")
    result = run_hifi_qc(BAM_FILE, "-o", output, "--sample-names", "MySample")
    assert result.returncode == 0, f"Failed: {result.stderr}"
    with open(output) as f:
        html = f.read()
    assert "MySample" in html
    print("PASS: custom_sample_name")


if __name__ == "__main__":
    test_basic_run()
    test_with_reference()
    test_sampling_mode()
    test_custom_sample_name()
    print("\nAll tests passed!")
