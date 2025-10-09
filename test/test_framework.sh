#!/bin/bash
#
# Simple bash testing framework for SPC_genome pipeline
#

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Arrays to store test results
declare -a FAILED_TESTS

#######################################
# Print colored output
#######################################
print_color() {
    local color=$1
    shift
    echo -e "${color}$@${NC}"
}

#######################################
# Assert that two values are equal
# Arguments:
#   $1: Expected value
#   $2: Actual value
#   $3: Message (optional)
#######################################
assert_equals() {
    local expected="$1"
    local actual="$2"
    local message="${3:-Values should be equal}"

    ((TESTS_RUN++))

    if [ "$expected" = "$actual" ]; then
        ((TESTS_PASSED++))
        print_color "$GREEN" "  ✓ PASS: $message"
        return 0
    else
        ((TESTS_FAILED++))
        FAILED_TESTS+=("$message")
        print_color "$RED" "  ✗ FAIL: $message"
        print_color "$RED" "    Expected: $expected"
        print_color "$RED" "    Actual:   $actual"
        return 1
    fi
}

#######################################
# Assert that a value is true (non-zero)
# Arguments:
#   $1: Value to test
#   $2: Message (optional)
#######################################
assert_true() {
    local value="$1"
    local message="${2:-Value should be true}"

    assert_equals "true" "$([[ $value ]] && echo true || echo false)" "$message"
}

#######################################
# Assert that a file exists
# Arguments:
#   $1: File path
#   $2: Message (optional)
#######################################
assert_file_exists() {
    local file="$1"
    local message="${2:-File $file should exist}"

    ((TESTS_RUN++))

    if [ -f "$file" ]; then
        ((TESTS_PASSED++))
        print_color "$GREEN" "  ✓ PASS: $message"
        return 0
    else
        ((TESTS_FAILED++))
        FAILED_TESTS+=("$message")
        print_color "$RED" "  ✗ FAIL: $message"
        return 1
    fi
}

#######################################
# Assert that a directory exists
# Arguments:
#   $1: Directory path
#   $2: Message (optional)
#######################################
assert_dir_exists() {
    local dir="$1"
    local message="${2:-Directory $dir should exist}"

    ((TESTS_RUN++))

    if [ -d "$dir" ]; then
        ((TESTS_PASSED++))
        print_color "$GREEN" "  ✓ PASS: $message"
        return 0
    else
        ((TESTS_FAILED++))
        FAILED_TESTS+=("$message")
        print_color "$RED" "  ✗ FAIL: $message"
        return 1
    fi
}

#######################################
# Assert that a command succeeds (exit code 0)
# Arguments:
#   $@: Command to run
#######################################
assert_success() {
    local message="Command should succeed: $*"

    ((TESTS_RUN++))

    if "$@" &>/dev/null; then
        ((TESTS_PASSED++))
        print_color "$GREEN" "  ✓ PASS: $message"
        return 0
    else
        ((TESTS_FAILED++))
        FAILED_TESTS+=("$message")
        print_color "$RED" "  ✗ FAIL: $message"
        return 1
    fi
}

#######################################
# Assert that a command fails (non-zero exit code)
# Arguments:
#   $@: Command to run
#######################################
assert_failure() {
    local message="Command should fail: $*"

    ((TESTS_RUN++))

    if ! "$@" &>/dev/null; then
        ((TESTS_PASSED++))
        print_color "$GREEN" "  ✓ PASS: $message"
        return 0
    else
        ((TESTS_FAILED++))
        FAILED_TESTS+=("$message")
        print_color "$RED" "  ✗ FAIL: $message"
        return 1
    fi
}

#######################################
# Print test summary
#######################################
print_summary() {
    echo ""
    echo "========================================"
    print_color "$YELLOW" "Test Summary"
    echo "========================================"
    echo "Total tests run:    $TESTS_RUN"
    print_color "$GREEN" "Passed:            $TESTS_PASSED"
    print_color "$RED" "Failed:            $TESTS_FAILED"

    if [ $TESTS_FAILED -gt 0 ]; then
        echo ""
        print_color "$RED" "Failed tests:"
        for test in "${FAILED_TESTS[@]}"; do
            print_color "$RED" "  - $test"
        done
    fi

    echo "========================================"

    # Return non-zero if any tests failed
    [ $TESTS_FAILED -eq 0 ]
}

#######################################
# Run all test functions in the current script
# Automatically detects functions starting with "test_"
#######################################
run_tests() {
    local test_file="${BASH_SOURCE[1]}"
    local test_name=$(basename "$test_file" .sh)

    print_color "$YELLOW" "Running: $test_name"
    echo ""

    # Find all functions starting with "test_"
    local test_functions=$(declare -F | awk '{print $3}' | grep '^test_')

    if [ -z "$test_functions" ]; then
        print_color "$YELLOW" "  No test functions found (functions should start with 'test_')"
        return 0
    fi

    # Run each test function
    for test_func in $test_functions; do
        echo "  Running: $test_func"
        $test_func
    done

    print_summary
}

#######################################
# Setup function - called before tests
# Override this in your test file
#######################################
setup() {
    :
}

#######################################
# Teardown function - called after tests
# Override this in your test file
#######################################
teardown() {
    :
}

# Export functions so they can be used in test files
export -f assert_equals
export -f assert_true
export -f assert_file_exists
export -f assert_dir_exists
export -f assert_success
export -f assert_failure
export -f print_summary
export -f print_color
export -f run_tests
