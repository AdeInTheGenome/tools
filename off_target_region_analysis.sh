#!/bin/bash
# Developer: Adetola Abdulkadir, SRS Bioinformatician II 

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  off_target_bcftools.sh -s <samples_dir> -b <bed_files_dir> -o <output_csv> [-d "1000,5000,10000,50000,1000000"] [--vd] [-q] [-h]

Required:
  -s   Path to samples directory (each sample has its own subdir)
  -b   Path to directory with .bed files
  -o   Output CSV filepath

Options:
  -d   Comma-separated distances to test (default: 1000,5000,10000,50000,1000000)
  --vd Run bcftools +variant-distance and add MedianNearestVarDist column
  -q   Quiet: suppress missing-VCF warnings
  -h   Help

Notes:
  - Requires: bcftools (and tabix), BED files are 0-based, half-open.
  - If using --vd, ensure the bcftools plugins are available. You can check with:
      bcftools plugin -l
    and set BCFTOOLS_PLUGINS if needed:
      export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
EOF
}

# Defaults
QUIET=0
WITH_VD=0
distances="1000,5000,10000,50000,1000000"

# Parse flags
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s) samples_dir="$2"; shift 2 ;;
    -b) bed_files_dir="$2"; shift 2 ;;
    -o) output_file="$2"; shift 2 ;;
    -d) distances="$2"; shift 2 ;;
    --vd) WITH_VD=1; shift ;;
    -q) QUIET=1; shift ;;
    -h) usage; exit 0 ;;
    *)  echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

: "${samples_dir:?Missing -s <samples_dir>}"
: "${bed_files_dir:?Missing -b <bed_files_dir>}"
: "${output_file:?Missing -o <output_csv>}"

command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found."; exit 4; }
command -v tabix >/dev/null 2>&1 || { echo "Error: tabix not found."; exit 4; }

[[ -d "$samples_dir" ]] || { echo "Error: samples_dir not found: $samples_dir"; exit 3; }
[[ -d "$bed_files_dir" ]] || { echo "Error: bed_files_dir not found: $bed_files_dir"; exit 3; }

shopt -s nullglob
shopt -s dotglob

# Collect .bed files
bed_list=("$bed_files_dir"/*.bed)
(( ${#bed_list[@]} > 0 )) || { echo "Error: No .bed files in $bed_files_dir"; exit 5; }

# Prepare output
if (( WITH_VD == 1 )); then
  echo "Sample,BED_File,Callset,Distance,Total_Mutations,PASS_Mutations,MedianNearestVarDist" > "$output_file"
else
  echo "Sample,BED_File,Callset,Distance,Total_Mutations,PASS_Mutations" > "$output_file"
fi

# Convert distance string to array
IFS=',' read -r -a DIST_ARR <<< "$distances"

found_any_sample=0

for sample_folder in "$samples_dir"/*/; do
  [[ -d "$sample_folder" ]] || continue
  found_any_sample=1
  sample_name=$(basename "$sample_folder")

  for bed_file in "${bed_list[@]}"; do
    bed_file_name=$(basename "$bed_file")

    for callset in sv small cnv; do
      vcf_file="${sample_folder%/}/${sample_name}_${callset}_off_target_edits.vcf.gz"

      if [[ ! -f "$vcf_file" ]]; then
        (( QUIET == 1 )) || echo "Warning: VCF not found: $vcf_file" >&2
        continue
      fi
      # Ensure indexed
      if [[ ! -f "${vcf_file}.tbi" && ! -f "${vcf_file}.csi" ]]; then
        tabix -p vcf "$vcf_file"
      fi

      for distance in "${DIST_ARR[@]}"; do
        # Create an expanded BED (±distance). Prevent negative starts.
        tmp_bed=$(mktemp)
        awk -v W="$distance" 'BEGIN{OFS="\t"} {s=$2-W; if(s<0)s=0; e=$3+W; print $1,s,e}' "$bed_file" > "$tmp_bed"

        # Total variants overlapping expanded regions
        total_count=$(bcftools view -R "$tmp_bed" "$vcf_file" -H | wc -l | awk '{print $1}')
        # PASS-only
        pass_count=$(bcftools view -R "$tmp_bed" -f PASS "$vcf_file" -H | wc -l | awk '{print $1}')

        if (( WITH_VD == 1 )); then
          # Annotate with +variant-distance, then compute median BasesToClosestVariant among the selected variants
          # If plugin missing, this will error; you can test availability with: bcftools plugin -l
          median_near=$(bcftools view -R "$tmp_bed" "$vcf_file" \
            | bcftools +variant-distance \
            | bcftools query -f '%INFO/BasesToClosestVariant\n' \
            | awk 'NF{a[NR]=$1} END{if(NR==0){print "NA"; exit}; asort(a); mid=int((NR+1)/2); if(NR%2){print a[mid]} else {print (a[mid]+a[mid+1])/2}}' \
            || echo "NA")
          echo "$sample_name,$bed_file_name,$callset,$distance,$total_count,$pass_count,$median_near" >> "$output_file"
        else
          echo "$sample_name,$bed_file_name,$callset,$distance,$total_count,$pass_count" >> "$output_file"
        fi

        rm -f "$tmp_bed"
      done
    done
  done
done

if (( found_any_sample == 0 )); then
  echo "Warning: No sample subdirectories found in: $samples_dir" >&2
fi

echo "bcftools-based off-target analysis completed. Results saved to $output_file."

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  off_target_analysis.sh -s <samples_dir> -b <bed_files_dir> -o <output_file> [-q]

Required:
  -s  Path to samples directory (each sample has its own folder)
  -b  Path to directory containing .bed files
  -o  Output CSV filepath

Options:
  -q  Quiet mode (suppress "missing VCF" warnings)
  -h  Show this help and exit

Example:
  off_target_analysis.sh -s adapted_vs_edited_samples -b bed_files -o off_target_report.csv
EOF
}

# Defaults
QUIET=0

# Parse flags
while getopts ":s:b:o:qh" opt; do
  case "$opt" in
    s) samples_dir="$OPTARG" ;;
    b) bed_files_dir="$OPTARG" ;;
    o) output_file="$OPTARG" ;;
    q) QUIET=1 ;;
    h) usage; exit 0 ;;
    \?) echo "Error: Invalid option -$OPTARG" >&2; usage; exit 2 ;;
    :)  echo "Error: Option -$OPTARG requires an argument." >&2; usage; exit 2 ;;
  esac
done

# Ensure required args
: "${samples_dir:?Missing -s <samples_dir>}"
: "${bed_files_dir:?Missing -b <bed_files_dir>}"
: "${output_file:?Missing -o <output_file>}"

# Basic validations
if [[ ! -d "$samples_dir" ]]; then
  echo "Error: samples_dir not found or not a directory: $samples_dir" >&2
  exit 3
fi
if [[ ! -d "$bed_files_dir" ]]; then
  echo "Error: bed_files_dir not found or not a directory: $bed_files_dir" >&2
  exit 3
fi
if ! command -v bedtools >/dev/null 2>&1; then
  echo "Error: bedtools not found in PATH." >&2
  exit 4
fi

# Safer globbing
shopt -s nullglob
shopt -s dotglob

# Collect .bed files; fail early if none
bed_list=("$bed_files_dir"/*.bed)
if (( ${#bed_list[@]} == 0 )); then
  echo "Error: No .bed files found in: $bed_files_dir" >&2
  exit 5
fi

# Initialize the output file with a header (truncate if exists)
echo "Sample,BED_File,Callset,Distance,Total_Mutations,PASS_Mutations" > "$output_file"

# Distance ranges
distances=(1000 5000 10000 50000 1000000)

# Looping over each sample folder (directories only)
found_any_sample=0
for sample_folder in "$samples_dir"/*/; do
  [[ -d "$sample_folder" ]] || continue
  found_any_sample=1
  sample_name=$(basename "$sample_folder")

  for bed_file in "${bed_list[@]}"; do
    bed_file_name=$(basename "$bed_file")

    # Looping over each variant callset
    for callset in sv small cnv; do
      vcf_file="${sample_folder%/}/${sample_name}_${callset}_off_target_edits.vcf.gz"

      # Confirming the VCF file exists
      if [[ ! -f "$vcf_file" ]]; then
        (( QUIET == 1 )) || echo "Warning: VCF file not found: $vcf_file" >&2
        continue
      fi

      # Looping through each specified distance
      for distance in "${distances[@]}"; do
        # Count all mutations within the specified distance
        total_count=$(bedtools window -a "$vcf_file" -b "$bed_file" -w "$distance" | wc -l)

        # Count only "PASS" mutations within the specified distance
        pass_count=$(bedtools window -a "$vcf_file" -b "$bed_file" -w "$distance" | awk '$7 == "PASS"' | wc -l)

        # Append results to the CSV file
        echo "$sample_name,$bed_file_name,$callset,$distance,$total_count,$pass_count" >> "$output_file"
      done
    done
  done
done

if (( found_any_sample == 0 )); then
  echo "Warning: No sample subdirectories found in: $samples_dir" >&2
fi

echo "Off-target analysis completed. Results saved to $output_file."
