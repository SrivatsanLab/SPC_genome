#!/bin/bash
#
# Example unit test demonstrating the test framework
#

# Source the test framework
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_framework.sh"

#######################################
# Test basic string equality
#######################################
test_string_equality() {
    local expected="hello"
    local actual="hello"
    assert_equals "$expected" "$actual" "Strings should match"
}

#######################################
# Test that config.yaml exists
#######################################
test_config_exists() {
    local config_file="$SCRIPT_DIR/../../config.yaml"
    assert_file_exists "$config_file" "config.yaml should exist in project root"
}

#######################################
# Test that scripts directory exists
#######################################
test_scripts_dir_exists() {
    local scripts_dir="$SCRIPT_DIR/../../scripts"
    assert_dir_exists "$scripts_dir" "scripts/ directory should exist"
}

#######################################
# Test that CapWGS_PP.sh is executable
#######################################
test_main_script_executable() {
    local main_script="$SCRIPT_DIR/../../CapWGS_PP.sh"
    assert_file_exists "$main_script" "CapWGS_PP.sh should exist"

    if [ -x "$main_script" ]; then
        print_color "$GREEN" "  ✓ PASS: CapWGS_PP.sh is executable"
        ((TESTS_PASSED++))
    else
        print_color "$RED" "  ✗ FAIL: CapWGS_PP.sh is not executable"
        ((TESTS_FAILED++))
    fi
    ((TESTS_RUN++))
}

# Run all tests
run_tests
