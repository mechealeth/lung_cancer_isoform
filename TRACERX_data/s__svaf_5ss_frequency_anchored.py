#!/usr/bin/env python3
import argparse
import os
import glob
import pysam


def iter_junction_bounds_1based(read):
    """
    Yield (left_1based, right_1based) for each spliced junction in a read.

    Definitions are purely genomic (independent of transcript strand):
      left  = last exonic base before the intron (lower coordinate side)
      right = first exonic base after the intron (higher coordinate side)

    For transcript-aware donor/acceptor:
      + strand: donor=left,  acceptor=right
      - strand: donor=right, acceptor=left
    """
    if read.is_unmapped or read.cigartuples is None:
        return

    ref = read.reference_start + 1  # 1-based reference cursor

    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X consume reference
            ref += length
        elif op == 2:        # D consumes reference
            ref += length
        elif op == 3:        # N (intron) consumes reference
            intron_start = ref                # 1-based first intronic base
            intron_end = ref + length - 1     # 1-based last intronic base
            left = intron_start - 1           # last exonic base before intron
            right = intron_end + 1            # first exonic base after intron
            yield left, right
            ref += length
        else:
            # I, S, H, P do not consume reference
            pass


def aligned_bp_overlap_with_exon(read, exon_start_1based, exon_end_1based):
    """
    Return number of aligned reference bases of 'read' that overlap
    [exon_start_1based, exon_end_1based] (1-based inclusive).

    Uses read.get_blocks(): list of (start0, end0) reference-aligned blocks,
    0-based half-open.
    """
    if read.is_unmapped:
        return 0

    exon_start0 = exon_start_1based - 1
    exon_end0 = exon_end_1based  # half-open end

    overlap = 0
    for b0, e0 in read.get_blocks():  # 0-based half-open
        s = max(b0, exon_start0)
        e = min(e0, exon_end0)
        if e > s:
            overlap += (e - s)
    return overlap


def count_svaf_frequency_anchored(
    bam_path,
    chrom,
    fixed_acceptor_1based,
    svaf_start_1based,
    svaf_end_1based,
    exon_start_1based,
    exon_end_1based,
    exon_anchor_bp=10,
    min_mapq=20,
    strand=None,          # "+", "-", or None
    count_mode="junction" # "junction" or "read"
):
    """
    Counts splicing INTO a fixed transcript 3' acceptor (fixed_acceptor_1based),
    using an exon-anchored definition:

      - read overlaps the given exon by >= exon_anchor_bp aligned bp
      - read contains splice(s) (CIGAR N)
      - (optional) read orientation matches the specified transcript strand

    Returns (J_total, J_svaf):
      - J_total: all qualifying junctions (or reads) that splice into fixed 3'SS
      - J_svaf: subset where the transcript donor (5'SS) lies within SVAF locus

    IMPORTANT (strand correctness):
      Junction bounds are computed in genomic coordinates (left/right).
      Then donor/acceptor are assigned based on transcript strand:
        + strand: donor=left,  acceptor=right
        - strand: donor=right, acceptor=left
      If strand=None, default assumes + strand donor/acceptor assignment.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Fetch reads overlapping the exon (anchor-based search space)
    exon_start0 = max(0, exon_start_1based - 1)
    exon_end0 = exon_end_1based

    J_total = 0
    J_svaf = 0

    # For count_mode="read", avoid counting the same read multiple times
    seen_total_reads = set()
    seen_svaf_reads = set()

    for read in bam.fetch(chrom, exon_start0, exon_end0):
        if read.is_unmapped:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.cigarstring is None or "N" not in read.cigarstring:
            continue

        # robustness: ignore secondary/supplementary/duplicates
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate:
            continue

        # Strand filter (orientation) if user specifies transcript strand
        # pysam: read.is_reverse True means read aligned to reverse strand
        if strand == "+" and read.is_reverse:
            continue
        if strand == "-" and (not read.is_reverse):
            continue

        # Require exon anchoring
        exon_olap = aligned_bp_overlap_with_exon(read, exon_start_1based, exon_end_1based)
        if exon_olap < exon_anchor_bp:
            continue

        hit_total = False
        hit_svaf = False

        for left, right in iter_junction_bounds_1based(read):
            # Assign transcript-aware donor/acceptor
            if strand == "-":
                donor = right
                acceptor = left
            else:
                # strand == "+" OR strand is None -> default (+) convention
                donor = left
                acceptor = right

            if acceptor == fixed_acceptor_1based:
                hit_total = True
                if svaf_start_1based <= donor <= svaf_end_1based:
                    hit_svaf = True

                if count_mode == "junction":
                    J_total += 1
                    if svaf_start_1based <= donor <= svaf_end_1based:
                        J_svaf += 1

        if count_mode == "read" and hit_total:
            # safer key for paired-end: distinguish read1/read2 if present
            # (for single-end, is_read1/is_read2 are False; that’s fine)
            qkey = (read.query_name, read.is_read1, read.is_read2)

            if qkey not in seen_total_reads:
                seen_total_reads.add(qkey)
                J_total += 1
            if hit_svaf and qkey not in seen_svaf_reads:
                seen_svaf_reads.add(qkey)
                J_svaf += 1

    bam.close()
    return J_total, J_svaf


def main():
    ap = argparse.ArgumentParser(
        description="Compute SVAF donor frequency into a fixed 3'SS using exon-anchored junction definition (strand-aware)."
    )
    ap.add_argument("--bam_glob", required=True)
    ap.add_argument("--chrom", required=True)

    ap.add_argument("--fixed_3ss", type=int, required=True,
                    help="Fixed acceptor coordinate (1-based) in TRANSCRIPT sense. "
                         "If --strand '-', the script will treat acceptor as the LOWER-coordinate junction boundary.")

    ap.add_argument("--svaf_start", type=int, required=True,
                    help="SVAF locus start (1-based, inclusive).")
    ap.add_argument("--svaf_end", type=int, required=True,
                    help="SVAF locus end (1-based, inclusive).")

    ap.add_argument("--exon_start", type=int, required=True,
                    help="Known exon start (1-based, inclusive).")
    ap.add_argument("--exon_end", type=int, required=True,
                    help="Known exon end (1-based, inclusive).")
    ap.add_argument("--exon_anchor_bp", type=int, default=10,
                    help="Require >= this many aligned bp overlapping the exon (default 10).")

    ap.add_argument("--min_mapq", type=int, default=20)
    ap.add_argument("--strand", choices=["+", "-"], default=None,
                    help="Transcript strand. Use '-' for negative-strand genes. "
                         "If unset, donor/acceptor assignment defaults to '+' convention.")
    ap.add_argument("--count_mode", choices=["junction", "read"], default="junction",
                    help="'junction' counts junction events; 'read' counts reads supporting the junction (default junction).")

    ap.add_argument("--out_tsv", required=True)

    args = ap.parse_args()

    bams = sorted(glob.glob(args.bam_glob))
    if not bams:
        raise SystemExit("No BAMs matched bam_glob")

    with open(args.out_tsv, "w") as out:
        out.write("sample\tbam\tJ_total\tJ_svaf\tf_svaf\n")
        for bam in bams:
            sample = os.path.basename(bam).replace(".Aligned.sortedByCoord.out.bam", "")
            J_total, J_svaf = count_svaf_frequency_anchored(
                bam_path=bam,
                chrom=args.chrom,
                fixed_acceptor_1based=args.fixed_3ss,
                svaf_start_1based=args.svaf_start,
                svaf_end_1based=args.svaf_end,
                exon_start_1based=args.exon_start,
                exon_end_1based=args.exon_end,
                exon_anchor_bp=args.exon_anchor_bp,
                min_mapq=args.min_mapq,
                strand=args.strand,
                count_mode=args.count_mode
            )
            f = (J_svaf / J_total) if J_total > 0 else 0.0
            out.write(f"{sample}\t{bam}\t{J_total}\t{J_svaf}\t{f:.6g}\n")


if __name__ == "__main__":
    main()
