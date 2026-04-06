#!/usr/bin/env python3
"""
Map per-chromosome est-sfs pval files back onto the VCF (preserving site order)
and produce:
 - {Output}_polarization.txt  (CHROM POS REF ALT REF_FREQ ALT_FREQ P_ANCESTRAL ANCESTRAL_ALLELE)
 - {Output}_ref_ancestral.bed (0-based half-open intervals for REF ancestral sites)

Assumptions:
 - est-sfs per-chromosome pval files contain no coordinates; they are in the same
   site order as the corresponding chromosome sites in the VCF.
 - pval file names include the chromosome name as the last '_',
   e.g. "PmuralisKimura-pval_CM014743.1.txt"  -> chrom = "CM014743.1"
"""

import pysam #type: ignore
import time
import os
import glob
from collections import deque, defaultdict

# ============== CONFIGURE ==============
VCF_PATH = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Pol_VCFs/Pre_polarization.vcf.gz"
PVAL_DIR = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Per_chr/Kimura"   # directory with per-chrom pval files
PVAL_GLOB = "*-pval*.txt"   # General pattern inside PVAL_DIR to find pval files
OUTPUT_PREFIX = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Per_chr/High_confidence"
OUTGROUPS = {"Psi_CAL_N84", "Pva_OU_E88", "SO08"}   # outgroup sample names to exclude when counting ingroup allele freqs
THRESHOLD = 0.95
# If True the script will raise a RuntimeError if ANY per-chromosome mismatch is found.
FAIL_ON_MISMATCH = True
# =======================================

def find_pval_files(pval_dir, pattern):
    files = sorted(glob.glob(os.path.join(pval_dir, pattern), recursive=True))
    return files

def pick_chrom_from_filename(fname, contigs):
    """Try to map filename -> contig using contigs list first, else fallback to last '_' token."""
    base = os.path.basename(fname)
    # try to find a contig name inside the filename (best)
    for c in contigs:
        if c in base:
            return c
    # fallback: remove final extension and take last '_' token
    if "." in base:
        base_noext = ".".join(base.split(".")[:-1])
    else:
        base_noext = base
    if "_" in base_noext:
        return base_noext.split("_")[-1]
    return base_noext

def load_pvals_by_chrom(files, contigs):
    """Return (pvals_by_chrom deque dict, pval_counts dict)"""
    pvals = defaultdict(deque)
    counts = {}
    if not files:
        return pvals, counts
    for f in files:
        chrom = pick_chrom_from_filename(f, contigs)
        n = 0
        with open(f, "rt") as fh:
            for line in fh:
                if not line.strip():
                    continue
                # skip est-sfs header lines that begin with '0'
                if line.startswith("0"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        p = float(parts[2])
                    except ValueError:
                        # skip unparsable lines (but continue reading)
                        continue
                    pvals[chrom].append(p)
                    n += 1
        counts[chrom] = counts.get(chrom, 0) + n
    return pvals, counts

def main():
    start_time = time.time()

    # verify paths
    if not os.path.exists(VCF_PATH):
        raise FileNotFoundError(f"VCF not found: {VCF_PATH}")
    if not os.path.isdir(PVAL_DIR):
        raise FileNotFoundError(f"PVAL_DIR not found or not a directory: {PVAL_DIR}")

    # open VCF to get contigs and samples
    vcf = pysam.VariantFile(VCF_PATH)
    contigs = list(vcf.header.contigs.keys())
    print(f"Loaded VCF: {VCF_PATH}")
    print(f"Contigs in VCF (first 10): {contigs[:10]}{'...' if len(contigs)>10 else ''}")

    # find pval files recursively
    files = find_pval_files(PVAL_DIR, PVAL_GLOB)
    if not files:
        # helpful debug listing
        print(f"No pval files found with pattern '{PVAL_GLOB}' under {PVAL_DIR}")
        raise FileNotFoundError(f"No pval files found under {PVAL_DIR} matching pattern {PVAL_GLOB}")

    print(f"Found {len(files)} pval files (example):")
    for f in files[:20]:
        print("  ", f)

    # load pvals mapped to chromosomes
    pvals_by_chrom, pval_counts = load_pvals_by_chrom(files, contigs)
    total_pvals = sum(len(dq) for dq in pvals_by_chrom.values())
    print(f"Loaded total pvals: {total_pvals}")
    print("Per-chrom pval counts (sample):")
    for k in sorted(pval_counts.keys()):
        print(f"  {k}: {pval_counts[k]}")

    # prepare samples and ingroup indices
    all_samples = list(vcf.header.samples)
    ingroup_samples = [s for s in all_samples if s not in OUTGROUPS]
    sample_indices = {s:i for i,s in enumerate(all_samples)}
    ingroup_indices = [sample_indices[s] for s in ingroup_samples]
    print(f"Ingroup samples (counted): {len(ingroup_indices)}  (outgroups excluded: {len(OUTGROUPS)})")

    # outputs
    txt_path = f"{OUTPUT_PREFIX}_polarization.txt"
    bed_path = f"{OUTPUT_PREFIX}_ref_ancestral.bed"
    os.makedirs(os.path.dirname(txt_path), exist_ok=True)
    txt_out = open(txt_path, "w")
    bed_out = open(bed_path, "w")
    txt_out.write("CHROM\tPOS\tREF\tALT\tREF_FREQ\tALT_FREQ\tP_ANCESTRAL\tANCESTRAL_ALLELE\n")

    # counters
    processed = skipped_non_snp = skipped_no_prob = 0
    high_conf_ref = 0
    used_pvals_by_chrom = defaultdict(int)
    total_vcf_snps_by_chrom = defaultdict(int)

    # iterate VCF and consume pvals per chromosome in order
    for record in vcf:
        # skip non-SNPs / multi-allelic
        if len(record.ref) != 1 or any(len(a) != 1 for a in record.alts):
            continue

        chrom = record.chrom
        # count this SNP in VCF per-chromosome (all SNPs irrespective of pval availability)
        total_vcf_snps_by_chrom[chrom] += 1

        if chrom not in pvals_by_chrom or not pvals_by_chrom[chrom]:
            # no pval available for this site/chrom
            skipped_no_prob += 1
            continue

        p_ancestral = pvals_by_chrom[chrom].popleft()
        used_pvals_by_chrom[chrom] += 1
        processed += 1

        # compute ingroup allele counts (0=REF, 1=ALT)
        ref_count = alt_count = 0
        for idx in ingroup_indices:
            sample_name = all_samples[idx]
            gt = record.samples[sample_name].get('GT', None)
            if gt is None or None in gt:
                continue
            for allele in gt:
                if allele == 0:
                    ref_count += 1
                elif allele == 1:
                    alt_count += 1

        total = ref_count + alt_count
        if total == 0:
            continue

        ref_freq = ref_count / total
        alt_freq = alt_count / total

        ancestral_allele = "unknown"
        if p_ancestral > THRESHOLD:
            if ref_freq >= 0.5:
                ancestral_allele = "REF"
                bed_out.write(f"{chrom}\t{record.pos-1}\t{record.pos}\n")
                high_conf_ref += 1
            else:
                ancestral_allele = "ALT"
        elif p_ancestral < (1-THRESHOLD):
            if ref_freq >= 0.5:
                ancestral_allele = "ALT"
            else:
                ancestral_allele = "REF"
                bed_out.write(f"{chrom}\t{record.pos-1}\t{record.pos}\n")
                high_conf_ref += 1

        if ancestral_allele != "unknown":
            txt_out.write(f"{chrom}\t{record.pos}\t{record.ref}\t{record.alts[0]}\t"
                          f"{ref_freq:.4f}\t{alt_freq:.4f}\t{p_ancestral:.6f}\t{ancestral_allele}\n")

    txt_out.close()
    bed_out.close()

    # Summary and strict checks
    print(f"\nProcessed (SNPs with pvals consumed): {processed}")
    print(f"Skipped sites with no pval available (VCF SNPs without pval): {skipped_no_prob}")
    print(f"High-confidence REF ancestral sites written to BED: {high_conf_ref}")
    print(f"TXT: {txt_path}")
    print(f"BED: {bed_path}")

    # Build union of chromosomes to report
    all_chroms = set(list(pvals_by_chrom.keys())) | set(list(used_pvals_by_chrom.keys())) | set(list(total_vcf_snps_by_chrom.keys()))
    print("\nPer-chromosome summary (vcf_snps / pvals_loaded / pvals_used / pvals_remaining):")
    mismatch_chroms = []
    leftover_chroms = []
    shortage_chroms = []
    for chrom in sorted(all_chroms):
        vcf_snps = total_vcf_snps_by_chrom.get(chrom, 0)
        loaded = pval_counts.get(chrom, 0)
        used = used_pvals_by_chrom.get(chrom, 0)
        remaining = len(pvals_by_chrom.get(chrom, ()))
        print(f"  {chrom}: {vcf_snps} / {loaded} / {used} / {remaining}")
        if remaining > 0:
            leftover_chroms.append(chrom)
        if vcf_snps != loaded:
            shortage_chroms.append(chrom)
        if vcf_snps > loaded:
            # too few pvals for VCF
            mismatch_chroms.append((chrom, "VCF_has_more", vcf_snps, loaded))
        elif loaded > vcf_snps:
            # more pvals than VCF SNPs
            mismatch_chroms.append((chrom, "PVALS_have_more", vcf_snps, loaded))

    # If any mismatch and FAIL_ON_MISMATCH, raise error with details
    if FAIL_ON_MISMATCH and mismatch_chroms:
        msg_lines = ["Per-chromosome mismatch detected (format: chrom, type, vcf_snps, pvals_loaded):"]
        for t in mismatch_chroms:
            msg_lines.append("  " + ", ".join(map(str, t)))
        msg_lines.append("Use the printed per-chromosome summary above to inspect counts.")
        raise RuntimeError("\n".join(msg_lines))

    if FAIL_ON_MISMATCH and leftover_chroms:
        raise RuntimeError(f"Pval files contain leftover (unused) pvals for chromosomes: {leftover_chroms}. See per-chrom summary above.")

    print("\nNo per-chromosome mismatches detected.")
    print(f"\nRuntime: {time.time() - start_time:.1f}s")

if __name__ == "__main__":
    main()