"""Generate a small synthetic HiFi BAM file for testing."""
import pysam
import random
import os

random.seed(42)

CHROMS = [("chr1", 10000), ("chr2", 8000)]
N_READS = 100
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def random_sequence(length):
    return "".join(random.choice("ACGT") for _ in range(length))


def create_test_bam():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    bam_path = os.path.join(OUTPUT_DIR, "test.bam")
    ref_path = os.path.join(OUTPUT_DIR, "test_ref.fa")

    # Create reference FASTA
    with open(ref_path, "w") as f:
        for name, length in CHROMS:
            seq = random_sequence(length)
            f.write(f">{name}\n{seq}\n")
    pysam.faidx(ref_path)

    # Create BAM
    header = pysam.AlignmentHeader.from_references(
        [name for name, _ in CHROMS],
        [length for _, length in CHROMS],
    )

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for i in range(N_READS):
            a = pysam.AlignedSegment()
            a.query_name = f"read_{i:04d}"
            chrom_idx = random.randint(0, len(CHROMS) - 1)
            chrom_name, chrom_len = CHROMS[chrom_idx]
            a.reference_id = chrom_idx

            # HiFi-like read length (1000-5000bp for our tiny genome)
            read_len = random.randint(1000, 5000)

            # Ensure read fits within chromosome
            max_pos = max(0, chrom_len - read_len - 50)
            a.reference_start = random.randint(0, max_pos)

            a.query_sequence = random_sequence(read_len)
            a.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(random.randint(30, 40) + 33) for _ in range(read_len))
            )

            # CIGAR: mostly matches with small indels
            cigar = []
            remaining = read_len
            while remaining > 0:
                op_len = min(random.randint(100, 1000), remaining)
                cigar.append((0, op_len))  # M
                remaining -= op_len
                if remaining > 5 and random.random() < 0.1:
                    ins_len = random.randint(1, 3)
                    if ins_len <= remaining:
                        cigar.append((1, ins_len))  # I
                        remaining -= ins_len
                if remaining > 5 and random.random() < 0.1:
                    cigar.append((2, random.randint(1, 3)))  # D

            a.cigar = cigar
            a.mapping_quality = 60
            a.flag = 0

            # HiFi tags
            a.set_tag("np", random.randint(3, 30), "i")
            a.set_tag("rq", round(random.uniform(0.99, 0.9999), 4), "f")

            bam.write(a)

    # Sort and index
    sorted_path = bam_path.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_path, bam_path)
    os.replace(sorted_path, bam_path)
    pysam.index(bam_path)

    print(f"Created: {bam_path}")
    print(f"Created: {ref_path}")


if __name__ == "__main__":
    create_test_bam()
