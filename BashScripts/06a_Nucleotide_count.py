#!/usr/bin/env python
"""
VCF Nucleotide Counter for Polarisation (est-sfs) - Chromosome-wise with Outgroup Haploidization

This script processes VCF files chromosome by chromosome, counting nucleotides per variant across specified sample groups.
For outgroups, it haploidizes homozygous sites and treats heterozygous sites as missing data.

Key features:
- Processes one chromosome at a time
- First column (ingroup) is separated by tab
- Outgroups are space-separated
- Output format: "A,C,G,T" counts
- Handles phased/unphased genotypes
- Skips non-SNPs and non-ACGT bases
- For outgroups: haploidizes homozygous sites, treats heterozygous as missing (0,0,0,0)
"""
import sys
import pysam #type:ignore
from collections import defaultdict
import os

def load_sample_groups(group_files):
    """Load sample names from group files - Order of the files matter. It must match the tree"""
    groups = []
    group_names = ["ingroup", "outgroup1", "outgroup2", "outgroup3"]
    for idx, filename in enumerate(group_files):
        try:
            with open(filename) as f:
                samples = [l.strip() for l in f if l.strip()]
            if not samples:
                sys.stderr.write(f"Error: Empty group file - {filename}\n")
                sys.exit(1)
            groups.append(samples)
            sys.stderr.write(f"Loaded {len(samples)} samples for {group_names[idx]}\n")
        except FileNotFoundError:
            sys.stderr.write(f"Error: Group file not found - {filename}\n")
            sys.exit(1)
    return groups

def get_chromosomes_with_data(vcf_path, groups):
    """Get ONLY chromosomes that actually have SNP data with samples"""
    sys.stderr.write("Scanning chromosomes for valid SNP data...\n")
    
    try:
        vcf = pysam.VariantFile(vcf_path)
        chromosomes_with_data = []
        
        all_chromosomes = list(vcf.header.contigs)
        sys.stderr.write(f"Found {len(all_chromosomes)} contigs in VCF header\n")
        
        for chrom in all_chromosomes:
            if any(pattern in chrom.lower() for pattern in ['scaffold', 'random', 'unplaced', 'unlocalized']):
                sys.stderr.write(f"Skipping scaffold/random: {chrom}\n")
                continue
                
            try:
                found_valid_snp = False
                record_count = 0
                
                for record in vcf.fetch(chrom):
                    record_count += 1
                    if record_count > 100 and found_valid_snp:
                        break
                    
                    # Check if this is a valid SNP
                    ref = record.ref.upper()
                    if len(ref) != 1 or ref not in 'ACGT':
                        continue
                    alts = tuple(alt.upper() for alt in record.alts or ())
                    if any(len(a) != 1 or a not in 'ACGT' for a in alts):
                        continue
                    
                    # Check if we have any sample data for this SNP
                    for gi, sample_list in enumerate(groups):
                        for sample in sample_list:
                            gt = record.samples[sample].get('GT')
                            if gt and not any(allele is None for allele in gt):
                                found_valid_snp = True
                                break
                        if found_valid_snp:
                            break
                    
                    if found_valid_snp:
                        break
                
                if found_valid_snp:
                    chromosomes_with_data.append(chrom)
                    sys.stderr.write(f"✓ {chrom}: Has valid SNP data ({record_count} variants scanned)\n")
                else:
                    sys.stderr.write(f"✗ {chrom}: No valid SNP data ({record_count} variants scanned)\n")
                    
            except ValueError:
                sys.stderr.write(f"✗ {chrom}: No variants found\n")
                continue
        
        vcf.close()
        return chromosomes_with_data
        
    except (IOError, ValueError) as e:
        sys.stderr.write(f"Error: Could not read VCF file {vcf_path}: {e}\n")
        sys.exit(1)

def process_chromosome(vcf_path, chromosome, groups):
    """Process a single chromosome and return nucleotide counts as lines
    
    Enhanced to handle multiallelic sites in outgroups while keeping ingroup biallelic.
    """
    try:
        vcf = pysam.VariantFile(vcf_path)
    except (IOError, ValueError) as e:
        sys.stderr.write(f"Error: Could not open VCF file {vcf_path}: {e}\n")
        return []

    output_lines = []
    variant_count = 0
    valid_variant_count = 0
    multiallelic_sites_count = 0
    
    try:
        for record in vcf.fetch(chromosome):
            variant_count += 1
            if variant_count % 10000 == 0:
                sys.stderr.write(f"  Processed {variant_count} variants on {chromosome}...\n")

            # Filter non-SNP and non-ACGT variants
            ref = record.ref.upper()
            if len(ref) != 1 or ref not in 'ACGT':
                continue
            alts = tuple(alt.upper() for alt in record.alts or ())
            if any(len(a) != 1 or a not in 'ACGT' for a in alts):
                continue

            # Check if this is multiallelic
            is_multiallelic = len(alts) > 1
            if is_multiallelic:
                multiallelic_sites_count += 1

            # Initialize nucleotide counters
            group_counts = [defaultdict(int) for _ in range(len(groups))]
            has_data = False

            # Process genotypes for each group
            for gi, sample_list in enumerate(groups):
                for sample in sample_list:
                    gt = record.samples[sample].get('GT')
                    
                    if not gt or any(allele is None for allele in gt):
                        continue
                        
                    # Ingroup processing (diploid) - MUST be biallelic for est-sfs
                    if gi == 0:
                        for allele in gt:
                            if allele == 0:
                                base = ref
                                group_counts[gi][base] += 1
                                has_data = True
                            elif allele == 1:
                                # For ingroup, only use first ALT (biallelic requirement)
                                if alts:
                                    base = alts[0]
                                    group_counts[gi][base] += 1
                                    has_data = True
                            # Skip alleles >1 in ingroup (est-sfs biallelic requirement)
                    
                    # Outgroup processing (haploidize with full multiallelic support)
                    else:
                        # Check if genotype is homozygous (all alleles the same)
                        if all(allele == gt[0] for allele in gt):
                            # Homozygous - count the allele
                            allele = gt[0]
                            if allele == 0:
                                base = ref
                                if base in 'ACGT':
                                    group_counts[gi][base] += 1
                                    has_data = True
                            elif allele > 0:
                                # For outgroups, use the actual ALT corresponding to the allele index
                                alt_index = allele - 1
                                if alt_index < len(alts):
                                    base = alts[alt_index]
                                    if base in 'ACGT':
                                        group_counts[gi][base] += 1
                                        has_data = True
                        # Heterozygous outgroups are treated as missing (do nothing)

            # Only include variant if it has data
            if has_data:
                ing = group_counts[0]
                ing_str = f"{ing['A']},{ing['C']},{ing['G']},{ing['T']}"
                out_strs = []
                for cnts in group_counts[1:]:
                    out_strs.append(f"{cnts['A']},{cnts['C']},{cnts['G']},{cnts['T']}")
                
                output_lines.append(f"{ing_str}\t{' '.join(out_strs)}\n")
                valid_variant_count += 1

    except ValueError as e:
        sys.stderr.write(f"Error processing chromosome {chromosome}: {e}\n")
        return []
    
    sys.stderr.write(f"  Chromosome {chromosome}: {valid_variant_count}/{variant_count} valid variants ({multiallelic_sites_count} multiallelic)\n")
    return output_lines

def main():
    if len(sys.argv) != 8:
        sys.stderr.write(
            "Usage: python count_nucleotides.py <input.vcf.gz> <chromosome|ALL> "
            "<output_dir> <ingroup.txt> <outgroup1.txt> <outgroup2.txt> <outgroup3.txt>\n"
        )
        sys.stderr.write(
            "Output format: <ingroup_counts>\t<outgroup1_counts> "
            "<outgroup2_counts> <outgroup3_counts>\n"
        )
        sys.stderr.write(
            "Features: Full multiallelic support for outgroups (homozygous only), biallelic ingroup for est-sfs\n"
        )
        sys.exit(1)

    vcf_path = sys.argv[1]
    chromosome_arg = sys.argv[2]
    output_dir = sys.argv[3]
    group_files = sys.argv[4:8]
    
    os.makedirs(output_dir, exist_ok=True)
    groups = load_sample_groups(group_files)
    
    if chromosome_arg.upper() == "ALL":
        chromosomes = get_chromosomes_with_data(vcf_path, groups)
        if not chromosomes:
            sys.stderr.write("ERROR: No chromosomes with valid SNP data found!\n")
            sys.exit(1)
        sys.stderr.write(f"\nProcessing {len(chromosomes)} chromosomes with data: {', '.join(chromosomes)}\n")
    else:
        all_chromosomes_with_data = get_chromosomes_with_data(vcf_path, groups)
        if chromosome_arg not in all_chromosomes_with_data:
            sys.stderr.write(f"ERROR: Chromosome {chromosome_arg} has no valid SNP data!\n")
            sys.exit(1)
        chromosomes = [chromosome_arg]
        sys.stderr.write(f"Processing specific chromosome: {chromosome_arg}\n")
    
    total_files_created = 0
    total_variants_written = 0
    
    for chromosome in chromosomes:
        sys.stderr.write(f"\nProcessing chromosome: {chromosome}\n")
        
        output_lines = process_chromosome(vcf_path, chromosome, groups)
        
        if output_lines:
            output_file = os.path.join(output_dir, f"nucleotide_counts_{chromosome}.tsv")
            with open(output_file, 'w') as out_f:
                out_f.writelines(output_lines)
            
            file_size = os.path.getsize(output_file)
            sys.stderr.write(f"✓ Created {output_file} ({len(output_lines)} variants, {file_size} bytes)\n")
            
            total_files_created += 1
            total_variants_written += len(output_lines)
        else:
            sys.stderr.write(f"✗ No data for {chromosome} - skipping file creation\n")
    
    sys.stderr.write(f"\n=== PROCESSING COMPLETE ===\n")
    sys.stderr.write(f"Total files created: {total_files_created}\n")
    sys.stderr.write(f"Total variants written: {total_variants_written}\n")
    
    if total_files_created == 0:
        sys.stderr.write("WARNING: No output files created!\n")
        sys.exit(1)

if __name__ == '__main__':
    main()