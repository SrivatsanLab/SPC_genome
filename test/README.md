# Testing Framework

This directory contains tests for the SPC_genome pipeline.

## Structure

```
test/
├── README.md                  # This file
├── test_runner.sh            # Main test runner script
├── unit/                     # Unit tests for individual scripts/functions
│   ├── test_barcode_parsing.sh
│   ├── test_chunk_creation.sh
│   └── ...
├── integration/              # Integration tests for full pipeline
│   ├── test_preprocessing_pipeline.sh
│   └── ...
└── fixtures/                 # Test data and expected outputs
    ├── sample_reads/
    └── expected_outputs/
```

## Running Tests

### Run All Tests
```bash
./test/test_runner.sh
```

### Run Specific Test Suite
```bash
./test/test_runner.sh unit          # Run only unit tests
./test/test_runner.sh integration   # Run only integration tests
```

### Run Individual Test
```bash
./test/unit/test_barcode_parsing.sh
```

## Writing Tests

Tests should follow this structure:

```bash
#!/bin/bash
source "$(dirname "$0")/../test_framework.sh"

test_function_name() {
    # Arrange
    input="test_value"
    expected="expected_output"

    # Act
    result=$(your_function "$input")

    # Assert
    assert_equals "$expected" "$result" "Description of what failed"
}

# Run tests
run_tests
```

## Test Categories

### Unit Tests
- Test individual functions/scripts in isolation
- Fast execution (< 1 minute)
- No SLURM job submission
- Use mock data from fixtures/

### Integration Tests
- Test full pipeline or major components
- May take longer (minutes to hours)
- May submit SLURM jobs
- Use small test datasets (< 1GB)

## Guidelines

1. **Keep tests fast**: Unit tests should run in seconds
2. **Use fixtures**: Store test data in fixtures/ directory
3. **Clean up**: Tests should clean up temporary files
4. **Be explicit**: Clear test names and assertion messages
5. **Independence**: Tests should not depend on each other
