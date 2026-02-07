#!/bin/bash

##########################################################################################################################
# Utility function to auto-detect known sites files for GATK BQSR
# This function looks for dbSNP, known indels, and Mills indels files in the reference bundle/ subdirectory
#
# Usage:
#   source scripts/utils/detect_known_sites.sh
#   detect_known_sites "/path/to/reference/dir"
#   # After calling, use these variables: DBSNP, KNOWN_INDELS, MILLS_INDELS, KNOWN_SITES_ARGS
#
# Returns:
#   Sets environment variables: DBSNP, KNOWN_INDELS, MILLS_INDELS, KNOWN_SITES_ARGS
#   Returns 0 if at least one known sites file was found, 1 if none were found
##########################################################################################################################

detect_known_sites() {
    local REFERENCE_DIR="$1"

    # Initialize variables as empty
    DBSNP=""
    KNOWN_INDELS=""
    MILLS_INDELS=""
    KNOWN_SITES_ARGS=""

    local found_count=0

    echo "========================================"
    echo "Detecting known sites files for BQSR..."
    echo "========================================"
    echo "Reference directory: ${REFERENCE_DIR}"
    echo ""

    # Look for files in bundle/ subdirectory
    local BUNDLE_DIR="${REFERENCE_DIR}/bundle"

    if [ ! -d "${BUNDLE_DIR}" ]; then
        echo "WARNING: bundle/ subdirectory not found at: ${BUNDLE_DIR}"
        echo "BQSR will be skipped. To use BQSR, ensure known sites files are available."
        echo ""
        return 1
    fi

    echo "Searching in: ${BUNDLE_DIR}"
    echo ""

    # Look for dbSNP (prefer newer versions)
    # Try dbsnp_146, then 144, then 138
    for version in 146 144 138; do
        local dbsnp_pattern="${BUNDLE_DIR}/dbsnp_${version}*.vcf.gz"
        local dbsnp_file=$(ls ${dbsnp_pattern} 2>/dev/null | head -n 1)

        if [ -n "$dbsnp_file" ] && [ -f "$dbsnp_file" ]; then
            DBSNP="$dbsnp_file"
            echo "✓ Found dbSNP: ${DBSNP}"
            found_count=$((found_count + 1))
            break
        fi
    done

    if [ -z "$DBSNP" ]; then
        echo "⚠ WARNING: dbSNP file not found (searched for dbsnp_146/144/138*.vcf.gz)"
    fi

    # Look for known indels
    local known_indels_pattern="${BUNDLE_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz"
    if [ -f "$known_indels_pattern" ]; then
        KNOWN_INDELS="$known_indels_pattern"
        echo "✓ Found known indels: ${KNOWN_INDELS}"
        found_count=$((found_count + 1))
    else
        echo "⚠ WARNING: Known indels file not found at: ${known_indels_pattern}"
    fi

    # Look for Mills indels
    local mills_pattern="${BUNDLE_DIR}/Mills_and_1000G_gold_standard.indels*.vcf.gz"
    local mills_file=$(ls ${mills_pattern} 2>/dev/null | head -n 1)

    if [ -n "$mills_file" ] && [ -f "$mills_file" ]; then
        MILLS_INDELS="$mills_file"
        echo "✓ Found Mills indels: ${MILLS_INDELS}"
        found_count=$((found_count + 1))
    else
        echo "⚠ WARNING: Mills indels file not found (searched for Mills_and_1000G_gold_standard.indels*.vcf.gz)"
    fi

    echo ""
    echo "Found ${found_count} known sites file(s)"
    echo ""

    # Build GATK arguments string
    if [ -n "$DBSNP" ]; then
        KNOWN_SITES_ARGS="${KNOWN_SITES_ARGS} --known-sites ${DBSNP}"
    fi
    if [ -n "$KNOWN_INDELS" ]; then
        KNOWN_SITES_ARGS="${KNOWN_SITES_ARGS} --known-sites ${KNOWN_INDELS}"
    fi
    if [ -n "$MILLS_INDELS" ]; then
        KNOWN_SITES_ARGS="${KNOWN_SITES_ARGS} --known-sites ${MILLS_INDELS}"
    fi

    # Export variables so they're available to calling script
    export DBSNP
    export KNOWN_INDELS
    export MILLS_INDELS
    export KNOWN_SITES_ARGS

    if [ "$found_count" -eq 0 ]; then
        echo "WARNING: No known sites files found. BQSR will be skipped."
        echo ""
        return 1
    fi

    echo "========================================"
    echo ""
    return 0
}
