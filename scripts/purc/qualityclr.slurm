#!/bin/bash
#SBATCH --job-name=analyze_clr
#SBATCH --partition=cpu-g2
#SBATCH --account=tribblelab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=0-9
#SBATCH --chdir=/mmfs1/gscratch/tribblelab/hemionitidoid_longreads
#SBATCH --export=ALL
#SBATCH -o log/analyze_clr_%A_%a.out
#SBATCH -e log/analyze_clr_%A_%a.err

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
PROJECT_ROOT="/mmfs1/gscratch/tribblelab/hemionitidoid_longreads"
BATCH_DIR_PARENT="${PROJECT_ROOT}/input/poolafastq"
OUTPUT_DIR="${PROJECT_ROOT}/output/clr_analysis"
SUMMARY_FILE="${OUTPUT_DIR}/summary_all_batches.txt"

# ------------------------------------------------------------
# Array job setup - processes batches 1-10
# ------------------------------------------------------------
BATCH_NUM=$(printf "%03d" "$((SLURM_ARRAY_TASK_ID + 1))")
BATCH_FOLDER="${BATCH_DIR_PARENT}/batch_${BATCH_NUM}"
INPUT_FASTQ="${BATCH_FOLDER}/reads.fastq"
OUTPUT_FILE="${OUTPUT_DIR}/batch_${BATCH_NUM}_analysis.txt"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Analyzing batch ${BATCH_NUM}"
echo "Input: ${INPUT_FASTQ}"
echo "Output: ${OUTPUT_FILE}"

# Check if input file exists
if [ ! -f "${INPUT_FASTQ}" ]; then
    echo "ERROR: Input file not found: ${INPUT_FASTQ}"
    exit 1
fi

# ------------------------------------------------------------
# Per-batch analysis
# ------------------------------------------------------------
{
    echo "=================================================="
    echo "CLR Read Analysis - Batch ${BATCH_NUM}"
    echo "Input file: ${INPUT_FASTQ}"
    echo "Analysis date: $(date)"
    echo "=================================================="
    echo ""

    # Count total reads
    TOTAL=$(grep -c "^@" "${INPUT_FASTQ}")
    echo "Total reads: $TOTAL"

    # Analyze read lengths and quality scores
    awk '
    BEGIN {
        short=0; mid=0; long=0;
        short_q=0; mid_q=0; long_q=0;
    }
    NR%4==2 {
        len=length($0);
        getline; getline;  # skip + line and get quality
        
        # Calculate average quality for this read
        q_sum=0;
        for(i=1; i<=length($0); i++) {
            q_sum += (index(" !\"#$%&'\''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", substr($0,i,1)) - 2);
        }
        avg_q = q_sum / length($0);
        
        # Categorize by length
        if (len < 800) {
            short++;
            short_q += avg_q;
        } else if (len <= 1300) {
            mid++;
            mid_q += avg_q;
        } else {
            long++;
            long_q += avg_q;
        }
    }
    END {
        print "\n=== Length Distribution ===";
        printf "< 800 bp:      %8d reads (%.1f%%)  Avg Q: %.1f\n", short, 100*short/(short+mid+long), (short>0 ? short_q/short : 0);
        printf "800-1300 bp:   %8d reads (%.1f%%)  Avg Q: %.1f\n", mid, 100*mid/(short+mid+long), (mid>0 ? mid_q/mid : 0);
        printf "> 1300 bp:     %8d reads (%.1f%%)  Avg Q: %.1f\n", long, 100*long/(short+mid+long), (long>0 ? long_q/long : 0);
        print "\n=== Quality Score Reference ===";
        print "Q0  (!) = 90% accuracy";
        print "Q10 (+) = 90% accuracy";
        print "Q20 (5) = 99% accuracy";
        print "Q30 (?) = 99.9% accuracy";
        print "Q40 (I) = 99.99% accuracy";
    }
    ' "${INPUT_FASTQ}"

    echo ""
    echo "=== Detailed Length Percentiles ==="
    awk 'NR%4==2 {print length($0)}' "${INPUT_FASTQ}" | sort -n | awk '
    BEGIN {n=0}
    {a[n++]=$1; sum+=$1}
    END {
        print "Min:    " a[0];
        print "5%:     " a[int(n*0.05)];
        print "25% Q1: " a[int(n*0.25)];
        print "50% Med:" a[int(n*0.50)];
        print "75% Q3: " a[int(n*0.75)];
        print "95%:    " a[int(n*0.95)];
        print "Max:    " a[n-1];
        print "Mean:   " int(sum/n);
    }
    '

    echo ""
    echo "=================================================="
    echo "Analysis complete for batch ${BATCH_NUM}"
    echo "=================================================="

} > "${OUTPUT_FILE}"

# Also print to stdout for the slurm log
cat "${OUTPUT_FILE}"

echo ""
echo "Results saved to: ${OUTPUT_FILE}"

# ------------------------------------------------------------
# Summary step - runs only on the last array task (index 9)
# Waits for all other batch output files to exist before proceeding
# ------------------------------------------------------------
if [ "${SLURM_ARRAY_TASK_ID}" -eq 9 ]; then
    echo ""
    echo "Last array task - waiting for all batch outputs before summarizing..."

    # Wait up to 10 minutes for all other batch files to appear
    MAX_WAIT=600
    WAITED=0
    while true; do
        MISSING=0
        for i in $(seq -w 1 10); do
            [ ! -f "${OUTPUT_DIR}/batch_${i}_analysis.txt" ] && MISSING=$((MISSING + 1))
        done
        [ "${MISSING}" -eq 0 ] && break
        if [ "${WAITED}" -ge "${MAX_WAIT}" ]; then
            echo "WARNING: Timed out waiting for ${MISSING} batch file(s). Summary may be incomplete."
            break
        fi
        sleep 10
        WAITED=$((WAITED + 10))
    done

    echo ""
    echo "=================================================="
    echo "CLR Read Analysis Summary - All Batches"
    echo "Date: $(date)"
    echo "=================================================="
    echo ""

    {
        echo "Batch | Total Reads | <800bp | 800-1300bp | >1300bp | Median | Mean"
        echo "------|-------------|--------|------------|---------|--------|------"

        for i in $(seq -w 1 10); do
            if [ -f "${OUTPUT_DIR}/batch_${i}_analysis.txt" ]; then
                TOTAL=$(grep "Total reads:"  "${OUTPUT_DIR}/batch_${i}_analysis.txt" | awk '{print $3}')
                SHORT=$(grep "< 800 bp:"     "${OUTPUT_DIR}/batch_${i}_analysis.txt" | awk '{print $3}')
                MID=$(grep "800-1300 bp:"    "${OUTPUT_DIR}/batch_${i}_analysis.txt" | awk '{print $3}')
                LONG=$(grep "> 1300 bp:"     "${OUTPUT_DIR}/batch_${i}_analysis.txt" | awk '{print $3}')
                MEDIAN=$(grep "50% Med:"     "${OUTPUT_DIR}/batch_${i}_analysis.txt" | awk '{print $3}')
                MEAN=$(grep "Mean:"          "${OUTPUT_DIR}/batch_${i}_analysis.txt" | awk '{print $2}')

                printf "%5s | %11s | %6s | %10s | %7s | %6s | %4s\n" \
                    "$i" "$TOTAL" "$SHORT" "$MID" "$LONG" "$MEDIAN" "$MEAN"
            fi
        done

        echo ""
        echo "Individual batch reports available in: ${OUTPUT_DIR}/batch_*_analysis.txt"

    } | tee "${SUMMARY_FILE}"

    echo ""
    echo "Summary saved to: ${SUMMARY_FILE}"
fi

