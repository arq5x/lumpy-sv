"""Regression tests for the BEDPE sorting helper."""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest


SCRIPT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "scripts", "bedpe_sort.py")


class BedpeSortTest(unittest.TestCase):
    def test_genome_order_does_not_drop_unlisted_contigs(self):
        rows = [
            "chrUn\t30\t31\tchrUn\t40\t41\tunknown-late",
            "chr2\t20\t21\tchr2\t30\t31\tknown-second",
            "chr1\t10\t11\tchr1\t15\t16\tknown-first",
            "chrUn\t5\t6\tchrUn\t8\t9\tunknown-early",
            "alt\t7\t8\talt\t9\t10\tother-unlisted",
        ]
        expected = [rows[2], rows[1], rows[4], rows[3], rows[0]]
        directory = tempfile.mkdtemp()
        try:
            bedpe_path = os.path.join(directory, "input.bedpe")
            genome_path = os.path.join(directory, "genome.txt")
            with open(bedpe_path, "w") as bedpe_file:
                bedpe_file.write("\n".join(rows) + "\n")
            with open(genome_path, "w") as genome_file:
                genome_file.write("chr1\t100\nchr2\t100\n")

            process = subprocess.Popen(
                [sys.executable, SCRIPT, "--bedpe_file", bedpe_path,
                 "--genome_file", genome_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = process.communicate()
        finally:
            shutil.rmtree(directory)

        self.assertEqual(process.returncode, 0, stderr)
        self.assertEqual(stdout.splitlines(), expected)


if __name__ == "__main__":
    unittest.main()
