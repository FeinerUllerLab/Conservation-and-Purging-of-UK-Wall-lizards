#!/bin/bash

############################################## SET UP ##############################################

# Path to the input directory:
indir=/home/feiner/Projects/UKwallies/RecombHotSpots/RawMaps

# Path to the output directory:
outdir=/home/feiner/Projects/UKwallies/RecombHotSpots/Step1

# Suffix of the input files (assuming the input file names follow the structure
# <chromosome>_<suffix>.rmap):
insuffix=VO01_n60_p2_b50_w100

# Path to the genome file (a tab-delimited file with the structure <chromName><TAB><chromSize>):
genome_file=/home/feiner/Projects/UKwallies/scripts/Chr_Length_GenBank.txt

# Window size in bp:
window_size=2000

# Step size in bp (for non-overlapping windows, set the step size equal to the window size):
step_size=1000

# Flank size in bp:
flank_size=40000

####################################################################################################

echo "Step 1: running (`date`)."

# Create the output directory if it does not already exist:
mkdir -p ${outdir}

# Create a suffix for the output files:
outsuffix=$(awk -v window_size=${window_size} -v step_size=${step_size} -v flank_size=${flank_size} 'BEGIN {OFS = ""; print "w", window_size / 1000, "kb", "_s", step_size / 1000, "kb", "_f", flank_size / 1000, "kb"}')

# Loop over the chromosomes:
while read -r chromosome chromosome_length; do

    echo "Processing chromosome ${chromosome} (${chromosome_length} bp)..."

    awk -v window_size=${window_size} -v step_size=${step_size} -v flank_size=${flank_size} -v chromosome=${chromosome} -v chromosome_length=${chromosome_length} '
    BEGIN {
        
        # Set the output field separator to tab:
        OFS = "\t"

        # Print a header:
        print "chromosome", "window_start", "window_end", "mean_r", "n_snps",
        "upstream_flank_start", "upstream_flank_end", "downstream_flank_start", "downstream_flank_end",
        "flanks_mean_r", "flanks_n_snps"

        # Determine the number of steps that can be taken before the next window goes beyond the end
        # of the chromosome:
        n_steps = int((chromosome_length - window_size) / step_size)

        # Add one step for any remaining stretch of chromosome:
        if ((chromosome_length - window_size) % step_size != 0) {
            n_steps++
        }

        # Create an array for each window variable with one element for each window (each window
        # being identified by an index going from 0 to n_steps):
        for (window_index = 0; window_index <= n_steps; window_index++) {

            # Window start position:
            window_start[window_index] = window_index * step_size

            # Window end position (cap at the chromosome length):
            if (window_start[window_index] + window_size <= chromosome_length) {
                window_end[window_index] = window_start[window_index] + window_size
            }
            else {
                window_end[window_index] = chromosome_length
            }

            # Numerator and denominator of the weighted mean recombination rate of the focal window:
            r_numerator[window_index] = 0
            r_denominator[window_index] = 0

            # Number of SNPs in the focal window:
            n_snps[window_index] = 0

            # Upstream flank start position (cap at zero):
            if (window_start[window_index] - flank_size >= 0) {
                upstream_flank_start[window_index] = window_start[window_index] - flank_size
            }
            else {
                upstream_flank_start[window_index] = 0
            }

            # Upstream flank end position:
            upstream_flank_end[window_index] = window_start[window_index]

            # Downstream flank start position:
            downstream_flank_start[window_index] = window_end[window_index]

            # Downstream flank end position (cap at the chromosome length):
            if (window_end[window_index] + flank_size <= chromosome_length) {
                downstream_flank_end[window_index] = window_end[window_index] + flank_size
            }
            else {
                downstream_flank_end[window_index] = chromosome_length
            }

            # Numerator and denominator of the weighted mean recombination rate of the upstream and
            # downstream flanks:
            upstream_flank_r_numerator[window_index] = 0
            upstream_flank_r_denominator[window_index] = 0
            downstream_flank_r_numerator[window_index] = 0
            downstream_flank_r_denominator[window_index] = 0

            # Number of SNPs in the upstream and downstream flanks:
            upstream_flank_n_snps[window_index] = 0
            downstream_flank_n_snps[window_index] = 0

        }

    }
    {
        # Get the start position, end position and recombination rate of each interval:
        interval_start = $1
        interval_end = $2
        interval_r = $3

        # Determine the number of steps required to get to the first window that overlaps with the
        # interval:
            # The expression int((interval_start - window_size) / step_size) + 1 yields the number
            # of steps required to get to the first window whose end position is greater than the
            # interval start position. The expression does not hold when the interval starts within
            # the first window (interval_start < window_size). In those cases, set the number of
            # steps to 0.
            # Derivation:
            # 1) window_end > interval_start
            # 2) window_start + window_size > interval_start
            # 3) window_index * step_size + window_size > interval_start
            # 4) window_index * step_size > interval_start - window_size
            # 5) window_index > (interval_start - window_size) / step_size
        if (interval_start >= window_size) {
            steps_to_first_overlap = int((interval_start - window_size) / step_size) + 1
        }
        else {
            steps_to_first_overlap = 0
        }

        # Determine the number of steps required to get to the last window that overlaps with the
        # interval:
            # The expression int(interval_end / step_size) yields the number of steps required to
            # get to the last window whose start position is less than or equal to the interval end
            # position. If the interval end position is equal to the window start position, there is
            # no overlap, because end positions are non-inclusive. In such cases, which occur when
            # the interval end position is exactly divisible by the step size, subtract one step. In
            # all cases, cap the number of steps at n_steps.
            # Derivation:
            # 1) window_start <= interval_end
            # 2) window_index * step_size <= interval_end
            # 3) window_index <= interval_end / step_size
        if (interval_end % step_size != 0) {
            steps_to_last_overlap = (int(interval_end / step_size) <= n_steps ? int(interval_end / step_size) : n_steps)
        }
        else {
            steps_to_last_overlap = ((interval_end / step_size) - 1 <= n_steps ? (interval_end / step_size) - 1 : n_steps)
        }

        # Loop over the windows that overlap with the interval:
        for (window_index = steps_to_first_overlap; window_index <= steps_to_last_overlap; window_index++) {

            # Define the start and end positions of the overlap between the interval and the window:
            overlap_start = (interval_start >= window_start[window_index] ? interval_start : window_start[window_index])
            overlap_end = (interval_end <= window_end[window_index] ? interval_end : window_end[window_index])

            # Add the contribution of the interval to the weighted mean recombination rate of the
            # window:
            r_numerator[window_index] += interval_r * (overlap_end - overlap_start)
            r_denominator[window_index] += (overlap_end - overlap_start)

            # If the interval starts within the window, add the SNP that defines the interval start
            # position to the window SNP count:
            if (interval_start >= window_start[window_index]) {
                n_snps[window_index]++
            }

        }

        # UPSTREAM FLANKS

        # Restrict the upstream flank calculations to intervals that contribute to at least one
        # upstream flank (i.e., exclude intervals whose start position is greater than or equal to
        # the upstream flank end position of the last window):
        if (interval_start < upstream_flank_end[n_steps]) {

            # Determine the number of steps required to get to the first window whose upstream flank
            # overlaps with the interval:
                # The expression int(interval_start / step_size) + 1 yields the number of steps
                # required to get to the first window whose upstream flank end position is greater
                # than the interval start position.
                # Derivation:
                # 1) upstream_flank_end > interval_start
                # 2) window_start > interval_start
                # 3) window_index * step_size > interval_start
                # 4) window_index > interval_start / step_size
            steps_to_first_upstream_flank_overlap = int(interval_start / step_size) + 1

            # Determine the number of steps required to get to the last window whose upstream flank
            # overlaps with the interval:
                # The expression int((interval_end + flank_size) / step_size) yields the number of
                # steps required to get to the last window whose upstream flank start position is
                # less than or equal to the interval end position. If the interval end position is
                # equal to the upstream flank start position, there is no overlap, because end
                # positions are non-inclusive. In such cases, which occur when the sum of the
                # interval end position and the flank size is exactly divisible by the step size,
                # subtract one step. In all cases, cap the number of steps at n_steps.
                # Derivation:
                # 1) upstream_flank_start <= interval_end
                # 2) window_start - flank_size <= interval_end
                # 3) window_index * step_size - flank_size <= interval_end
                # 4) window_index * step_size <= interval_end + flank_size
                # 5) window_index <= (interval_end + flank_size) / step_size
            if ((interval_end + flank_size) % step_size != 0) {
                steps_to_last_upstream_flank_overlap = (int((interval_end + flank_size) / step_size) <= n_steps ? int((interval_end + flank_size) / step_size) : n_steps)
            }
            else {
                steps_to_last_upstream_flank_overlap = (int((interval_end + flank_size) / step_size) - 1 <= n_steps ? int((interval_end + flank_size) / step_size) - 1 : n_steps)
            }

            # Loop over the windows whose upstream flanks overlap with the interval:
            for (window_index = steps_to_first_upstream_flank_overlap; window_index <= steps_to_last_upstream_flank_overlap; window_index++) {

                # Define the start and end positions of the overlap between the interval and the
                # upstream flank:
                overlap_start = (interval_start >= upstream_flank_start[window_index] ? interval_start : upstream_flank_start[window_index])
                overlap_end = (interval_end <= upstream_flank_end[window_index] ? interval_end : upstream_flank_end[window_index])

                # Add the contribution of the interval to the weighted mean recombination rate of
                # the upstream flank:
                upstream_flank_r_numerator[window_index] += interval_r * (overlap_end - overlap_start)
                upstream_flank_r_denominator[window_index] += (overlap_end - overlap_start)

                # If the interval starts within the upstream flank, add the SNP that defines the
                # interval start position to the upstream flank SNP count:
                if (interval_start >= upstream_flank_start[window_index]) {
                    upstream_flank_n_snps[window_index]++
                }

            }

        }
    
        # DOWNSTREAM FLANKS

        # Restrict the downstream flank calculations to intervals that contribute to at least one
        # downstream flank (i.e., exclude intervals whose end position is less than or equal to the
        # downstream flank start position of the first window):
        if (interval_end > downstream_flank_start[0]) {
            
            # Determine the number of steps required to get to the first window whose downstream
            # flank overlaps with the interval:
                # The expression int((interval_start - window_size - flank_size) / step_size) + 1
                # yields the number of steps required to get to the first window whose downstream
                # flank end position is greater than the interval start position. The expression
                # does not hold when the interval starts within the downstream flank of the first
                # window (interval start < window size + flank size). In those cases, set the number
                # of steps to 0.
                # Derivation:
                # 1) downstream_flank_end > interval_start
                # 2) window_end + flank_size > interval_start
                # 3) window_start + window_size + flank_size > interval_start
                # 4) window_index * step_size + window_size + flank_size > interval_start
                # 5) window_index * step_size > interval_start - window_size - flank_size
                # 6) window_index > (interval_start - window_size - flank_size) / step_size
            if (interval_start >= window_size + flank_size) {
                steps_to_first_downstream_flank_overlap = int((interval_start - window_size - flank_size) / step_size) + 1
            }
            else {
                steps_to_first_downstream_flank_overlap = 0
            }

            # Determine the number of steps required to get to the last window whose downstream
            # flank overlaps with the interval:
                # The expression int((interval_end - window_size) / step_size) yields the number of
                # steps required to get to the last window whose downstream flank start position is
                # less than or equal to the interval end position. If the interval end position is
                # equal to the downstream flank start position, there is no overlap, because end
                # positions are non-inclusive. In such cases, which occur when the interval end
                # position minus the window size is exactly divisible by the step size, subtract one
                # step. In all cases, cap the number of steps at n_steps.
                # Derivation:
                # 1) downstream_flank_start <= interval_end
                # 2) window_end <= interval_end
                # 3) window_start + window_size <= interval_end
                # 4) window_index * step_size + window_size <= interval_end
                # 5) window_index * step_size <= interval_end - window_size
                # 6) window_index <= (interval_end - window_size) / step_size
            if ((interval_end - window_size) % step_size != 0) {
                steps_to_last_downstream_flank_overlap = (int((interval_end - window_size) / step_size) <= n_steps ? int((interval_end - window_size) / step_size) : n_steps)
            }
            else {
                steps_to_last_downstream_flank_overlap = (int((interval_end - window_size) / step_size) - 1 <= n_steps ? int((interval_end - window_size) / step_size) - 1 : n_steps)
            }

            # Loop over the windows whose downstream flanks overlap with the interval:
            for (window_index = steps_to_first_downstream_flank_overlap; window_index <= steps_to_last_downstream_flank_overlap; window_index++) {

                # Define the start and end positions of the overlap between the interval and the
                # downstream flank:
                overlap_start = (interval_start >= downstream_flank_start[window_index] ? interval_start : downstream_flank_start[window_index])
                overlap_end = (interval_end <= downstream_flank_end[window_index] ? interval_end : downstream_flank_end[window_index])

                # Add the contribution of the interval to the weighted mean recombination rate of
                # the downstream flank:
                downstream_flank_r_numerator[window_index] += interval_r * (overlap_end - overlap_start)
                downstream_flank_r_denominator[window_index] += (overlap_end - overlap_start)

                # If the interval starts within the downstream flank, add the SNP that defines the
                # interval start position to the downstream flank SNP count:
                if (interval_start >= downstream_flank_start[window_index]) {
                    downstream_flank_n_snps[window_index]++
                }

            }

        }
    }
    END {

        # Loop over all the windows:
        for (window_index = 0; window_index <= n_steps; window_index++) {

            # Calculate the weighted mean recombination rate of the window:
            mean_r[window_index] = (r_denominator[window_index] > 0 ? r_numerator[window_index] / r_denominator[window_index] : "NA")

            # Calculate the weighted mean recombination rate of the upstream flank:
            #upstream_flank_mean_r[window_index] = (upstream_flank_r_denominator[window_index] > 0 ? upstream_flank_r_numerator[window_index] / upstream_flank_r_denominator[window_index] : "NA")

            # Calculate the weighted mean recombination rate of the downstream flank:
            #downstream_flank_mean_r[window_index] = (downstream_flank_r_denominator[window_index] > 0 ? downstream_flank_r_numerator[window_index] / downstream_flank_r_denominator[window_index] : "NA")

            # Calculate the weighted mean recombination rate of both flanks combined:
            flanks_r_numerator[window_index] = upstream_flank_r_numerator[window_index] + downstream_flank_r_numerator[window_index]
            flanks_r_denominator[window_index] = upstream_flank_r_denominator[window_index] + downstream_flank_r_denominator[window_index]
            flanks_mean_r[window_index] = (flanks_r_denominator[window_index] > 0 ? flanks_r_numerator[window_index] / flanks_r_denominator[window_index] : "NA")

            # Count the number of SNPs in both flanks combined:
            flanks_n_snps[window_index] = upstream_flank_n_snps[window_index] + downstream_flank_n_snps[window_index]

            # Print the window information:
            print chromosome, window_start[window_index], window_end[window_index], mean_r[window_index], n_snps[window_index],
            upstream_flank_start[window_index], upstream_flank_end[window_index], downstream_flank_start[window_index], downstream_flank_end[window_index],
            flanks_mean_r[window_index], flanks_n_snps[window_index]

        }

    }
    ' ${indir}/${chromosome}_${insuffix}.rmap > ${outdir}/${chromosome}_${outsuffix}.txt

done < ${genome_file}

echo "Step 1: done (`date`)."
