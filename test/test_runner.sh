#!/bin/bash
#
# Test runner for SPC_genome pipeline
# Runs all tests or specific test suites
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Test categories to run
RUN_UNIT=false
RUN_INTEGRATION=false

# Parse arguments
if [ $# -eq 0 ]; then
    # Run all tests by default
    RUN_UNIT=true
    RUN_INTEGRATION=true
else
    while [ $# -gt 0 ]; do
        case "$1" in
            unit)
                RUN_UNIT=true
                ;;
            integration)
                RUN_INTEGRATION=true
                ;;
            -h|--help)
                echo "Usage: $0 [unit|integration]"
                echo ""
                echo "Run tests for the SPC_genome pipeline"
                echo ""
                echo "Arguments:"
                echo "  unit          Run only unit tests"
                echo "  integration   Run only integration tests"
                echo "  (none)        Run all tests"
                echo ""
                echo "Examples:"
                echo "  $0              # Run all tests"
                echo "  $0 unit         # Run only unit tests"
                echo "  $0 integration  # Run only integration tests"
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
        shift
    done
fi

# Initialize counters
TOTAL_TESTS=0
TOTAL_PASSED=0
TOTAL_FAILED=0

echo -e "${YELLOW}================================${NC}"
echo -e "${YELLOW}SPC Genome Pipeline Test Suite${NC}"
echo -e "${YELLOW}================================${NC}"
echo ""

#######################################
# Run tests in a directory
#######################################
run_test_suite() {
    local suite_name=$1
    local test_dir=$2

    if [ ! -d "$test_dir" ]; then
        echo -e "${YELLOW}No $suite_name tests found in $test_dir${NC}"
        return 0
    fi

    local test_files=$(find "$test_dir" -name "test_*.sh" -type f 2>/dev/null || true)

    if [ -z "$test_files" ]; then
        echo -e "${YELLOW}No $suite_name test files found${NC}"
        echo ""
        return 0
    fi

    echo -e "${YELLOW}Running $suite_name tests...${NC}"
    echo ""

    local suite_passed=0
    local suite_failed=0

    for test_file in $test_files; do
        if [ -x "$test_file" ]; then
            if bash "$test_file"; then
                ((suite_passed++))
            else
                ((suite_failed++))
            fi
        else
            echo -e "${YELLOW}Warning: $test_file is not executable, skipping${NC}"
        fi
        echo ""
    done

    TOTAL_PASSED=$((TOTAL_PASSED + suite_passed))
    TOTAL_FAILED=$((TOTAL_FAILED + suite_failed))
    TOTAL_TESTS=$((TOTAL_TESTS + suite_passed + suite_failed))

    echo -e "${YELLOW}$suite_name summary: ${GREEN}$suite_passed passed${NC}, ${RED}$suite_failed failed${NC}"
    echo ""
}

# Run unit tests
if [ "$RUN_UNIT" = true ]; then
    run_test_suite "Unit" "$SCRIPT_DIR/unit"
fi

# Run integration tests
if [ "$RUN_INTEGRATION" = true ]; then
    run_test_suite "Integration" "$SCRIPT_DIR/integration"
fi

# Print final summary
echo -e "${YELLOW}================================${NC}"
echo -e "${YELLOW}Final Summary${NC}"
echo -e "${YELLOW}================================${NC}"
echo "Total tests:     $TOTAL_TESTS"
echo -e "${GREEN}Passed:         $TOTAL_PASSED${NC}"
echo -e "${RED}Failed:         $TOTAL_FAILED${NC}"
echo -e "${YELLOW}================================${NC}"

# Exit with appropriate code
if [ $TOTAL_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed${NC}"
    exit 1
fi
