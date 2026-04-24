#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./check_positions_in_vcfs.sh /path/to/vcf_dir positions.txt > results.tsv
#
# positions.txt:  CHROM <tab> POS

if [ $# -ne 2 ]; then
    echo "Usage: $0 /path/to/vcf_dir positions.txt" >&2
    exit 1
fi

VCF_DIR="$1"
POS_FILE="$2"

if [ ! -d "$VCF_DIR" ]; then
    echo "ERROR: VCF directory '$VCF_DIR' does not exist." >&2
    exit 1
fi

if [ ! -f "$POS_FILE" ]; then
    echo "ERROR: Position file '$POS_FILE' does not exist." >&2
    exit 1
fi

# Gather VCF files manually (no mapfile)
VCF_FILES=""
for f in "$VCF_DIR"/*.vcf.gz; do
    # Handle no matches
    [ -e "$f" ] || continue
    VCF_FILES="$VCF_FILES $f"
done

if [ -z "$VCF_FILES" ]; then
    echo "ERROR: No .vcf.gz files found in '$VCF_DIR'." >&2
    exit 1
fi

# Check for bcftools
if command -v bcftools >/dev/null 2>&1; then
    HAVE_BCFTOOLS=1
else
    HAVE_BCFTOOLS=0
    echo "WARNING: bcftools not found. Using slower zgrep/awk search." >&2
fi

echo -e "CHROM\tPOS\tFOUND\tVCF_FILES"

# Read the input file line-by-line
while IFS=$'\t' read -r CHROM POS _; do
    # Skip empty or commented lines
    [ -z "$CHROM" ] && continue
    case "$CHROM" in \#*) continue ;; esac

    FOUND_VCFS=""

    for VCF in $VCF_FILES; do
        if [ $HAVE_BCFTOOLS -eq 1 ]; then
            if bcftools view -H -r "${CHROM}:${POS}-${POS}" "$VCF" 2>/dev/null \
                | grep -q .; then
                FOUND_VCFS="${FOUND_VCFS},$(basename "$VCF")"
            fi
        else
            if zgrep -v '^#' "$VCF" 2>/dev/null \
                | awk -v C="$CHROM" -v P="$POS" '$1==C && $2==P {found=1} END{exit !found}'; then
                FOUND_VCFS="${FOUND_VCFS},$(basename "$VCF")"
            fi
        fi
    done

    if [ -n "$FOUND_VCFS" ]; then
        # Remove leading comma
        FOUND_VCFS="${FOUND_VCFS#,}"
        printf "%s\t%s\tyes\t%s\n" "$CHROM" "$POS" "$FOUND_VCFS"
    else
        printf "%s\t%s\tno\t-\n" "$CHROM" "$POS"
    fi

done < "$POS_FILE"
