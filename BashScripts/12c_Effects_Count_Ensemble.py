#!/usr/bin/env python3
"""
Objective: Count SnpEff impact categories from a VCF file. 
"""

import gzip
import collections

def count_impacts(vcf_path, output_path):
    impact_counts = collections.Counter()
    
    with gzip.open(vcf_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue  # Skip header lines
                
            fields = line.strip().split('\t')
            info_field = fields[7]  # INFO column
            
            # Extract IMPACT from ANN field
            for info_item in info_field.split(';'):
                if info_item.startswith('ANN='):
                    ann_value = info_item[4:]  # Remove 'ANN='
                    # Get first annotation's impact
                    first_ann = ann_value.split(',')[0]
                    impact = first_ann.split('|')[2]  # IMPACT is the 3rd field in ANN
                    impact_counts[impact] += 1
                    break  # Only count each site once
    
    # Totals
    total_all = sum(impact_counts.values())
    total_no_mod = sum(count for imp, count in impact_counts.items() if imp != "MODIFIER")
    
    # Write results
    with open(output_path, 'w') as out:
        out.write("IMPACT\tCOUNT\tPERCENT_TOTAL\tPERCENT_NO_MODIFIER\n")
        for impact, count in impact_counts.most_common():
            pct_total = (count / total_all * 100) if total_all > 0 else 0
            pct_no_mod = (count / total_no_mod * 100) if (impact != "MODIFIER" and total_no_mod > 0) else "-"
            out.write(f"{impact}\t{count}\t{pct_total:.2f}\t{pct_no_mod}\n")
    
    print(f"Counts saved to: {output_path}")
    for impact, count in impact_counts.most_common():
        pct_total = (count / total_all * 100) if total_all > 0 else 0
        pct_no_mod = (count / total_no_mod * 100) if (impact != "MODIFIER" and total_no_mod > 0) else "-"
        print(f"{impact}: {count} ({pct_total:.2f}% total, {pct_no_mod}% no-modifier)")

if __name__ == "__main__":
    vcf_path = "/home/feiner/Projects/UKwallies/Purging/All_Polarised_annotated_Ensemble.vcf.gz"
    output_path = "/home/feiner/Projects/UKwallies/Purging/Effect_Counts_Polarised_Final_Ensemble.tsv"
    count_impacts(vcf_path, output_path)
