package utils

import (
	"fmt"
	"numerical-methods/linalg"
	"os"
	"strconv"
	"strings"
)

func ParseRow(line string, n int) (linalg.Vector, error) {
	fields := strings.Fields(line)
	if len(fields) != n {
		return nil, fmt.Errorf("expected %d numbers, got %d", n, len(fields))
	}

	row := make(linalg.Vector, n)
	for i, field := range fields {
		value, err := strconv.ParseFloat(field, 64)
		if err != nil {
			return nil, fmt.Errorf("%q is not a valid number", field)
		}
		row[i] = value
	}

	return row, nil
}

func ParseEpsilon(inputPath string) (float64, error) {
	data, err := os.ReadFile(inputPath)
	if err != nil {
		return 0, fmt.Errorf("read input file: %w", err)
	}

	fields := strings.Fields(string(data))
	if len(fields) != 1 {
		return 0, fmt.Errorf("input file must contain exactly one epsilon value")
	}

	epsilon, err := strconv.ParseFloat(fields[0], 64)
	if err != nil {
		return 0, fmt.Errorf("parse epsilon: %w", err)
	}
	if epsilon <= 0 {
		return 0, fmt.Errorf("epsilon must be positive")
	}

	return epsilon, nil
}
