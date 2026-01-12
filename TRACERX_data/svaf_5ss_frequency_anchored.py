###Question :Among transcripts that splice into the fixed gene 3′SS, what fraction use any donor within the SVAF locus? ###

#!/usr/bin/env python3
import argparse
import os
import glob
import pysam

def iter_junctions_1based(read):
    """
    Yield (donor_1based, acceptor_1based) for each spliced junction in a read.
    donor = last exonic base before intron
    acceptor = last intronic base (boundary at intron end)
    """
    if read.is_unmapped or read.cigartuples is None:
        return
    ref = read.reference_start + 1  # 1-based cursor

    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X consume reference
            ref += length
        elif op == 2:        # D consumes reference
            ref += length
        elif op == 3:        # N (intron) consumes reference
            intron_start = ref                # 1-based first intronic base
            intron_end = ref + length - 1     # 1-based last intronic base
            donor = intron_start - 1
            acceptor = intron_end
            yield donor, acceptor
            ref += length
        else:
            # I, S, H, P do not consume reference
            pass

def aligned_bp_overlap_with_exon(read, exon_start_1based, exon_end_1based):
    """
    Return number of aligned reference bases of 'read' that overlap [exon_start, exon_end] (1-based inclusive).
    Uses read.get_blocks(): list of (start0, end0) reference-aligned blocks, 0-based half-open.
    """
    if read.is_unmapped:
        return 0
    exon_start0 = exon_start_1based - 1
    exon_end0 = exon_end_1based       # half-open end

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
    strand=None,
    count_mode="junction"  # "junction" or "read"
):
    """
    Counts junctions into fixed_acceptor_1based with additional requirement:
      read has >= exon_anchor_bp aligned bp overlapping [exon_start_1based, exon_end_1based].

    Returns (J_total, J_svaf).
    - J_total: all qualifying junctions (or reads) into fixed 3'SS
    - J_svaf: subset with donor within SVAF locus
    """

    bam = pysam.AlignmentFile(bam_path, "rb")

    # Fetch reads that overlap the exon (this is the anchor-based "search space")
    exon_start0 = max(0, exon_start_1based - 1)
    exon_end0 = exon_end_1based

    J_total = 0
    J_svaf = 0

    # For count_mode="read", avoid counting the same read twice
    # (if it has multiple junctions hitting the acceptor)
    seen_total_reads = set()
    seen_svaf_reads = set()

    for read in bam.fetch(chrom, exon_start0, exon_end0):
        if read.is_unmapped:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.cigarstring is None or "N" not in read.cigarstring:
            continue

        if strand == "+" and read.is_reverse:
            continue
        if strand == "-" and (not read.is_reverse):
            continue

        # Require exon anchoring
        exon_olap = aligned_bp_overlap_with_exon(read, exon_start_1based, exon_end_1based)
        if exon_olap < exon_anchor_bp:
            continue

        # Now examine junction(s) in this read
        hit_total = False
        hit_svaf = False
        for donor, acceptor in iter_junctions_1based(read):
            if acceptor == fixed_acceptor_1based:
                hit_total = True
                if svaf_start_1based <= donor <= svaf_end_1based:
                    hit_svaf = True
                if count_mode == "junction":
                    J_total += 1
                    if svaf_start_1based <= donor <= svaf_end_1based:
                        J_svaf += 1

        if count_mode == "read" and hit_total:
            qn = read.query_name
            if qn not in seen_total_reads:
                seen_total_reads.add(qn)
                J_total += 1
            if hit_svaf and qn not in seen_svaf_reads:
                seen_svaf_reads.add(qn)
                J_svaf += 1

    bam.close()
    return J_total, J_svaf

def main():
    ap = argparse.ArgumentParser(
        description="Compute SVAF donor frequency into a fixed 3'SS using exon-anchored junction definition."
    )
    ap.add_argument("--bam_glob", required=True)
    ap.add_argument("--chrom", required=True)

    ap.add_argument("--fixed_3ss", type=int, required=True,
                    help="Fixed acceptor coordinate (1-based).")

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
                    help="Optional: restrict to reads consistent with strand. Leave unset if unsure.")
    ap.add_argument("--count_mode", choices=["junction", "read"], default="junction",
                    help="'junction' counts junction events; 'read' counts reads that support the junction (default junction).")

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

