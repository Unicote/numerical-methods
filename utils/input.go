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
	rows, err := ReadRows(inputPath, 1)
	if err != nil {
		return 0, fmt.Errorf("parse epsilon: %w", err)
	}

	epsilon := rows[0][0]
	if epsilon <= 0 {
		return 0, fmt.Errorf("epsilon must be positive")
	}

	return epsilon, nil
}

func ReadLines(inputPath string) ([]string, error) {
	data, err := os.ReadFile(inputPath)
	if err != nil {
		return nil, fmt.Errorf("read input file: %w", err)
	}

	lines := make([]string, 0)
	for line := range strings.SplitSeq(string(data), "\n") {
		line = strings.TrimSpace(line)
		if line != "" {
			lines = append(lines, line)
		}
	}
	return lines, nil
}

func ReadRows(inputPath string, sizes ...int) ([]linalg.Vector, error) {
	lines, err := ReadLines(inputPath)
	if err != nil {
		return nil, err
	}
	if len(lines) != len(sizes) {
		return nil, fmt.Errorf("input file must contain %d non-empty lines, got %d", len(sizes), len(lines))
	}

	rows := make([]linalg.Vector, len(lines))
	for i, line := range lines {
		row, err := ParseRow(line, sizes[i])
		if err != nil {
			return nil, fmt.Errorf("line %d of input file: %w", i+1, err)
		}
		rows[i] = row
	}
	return rows, nil
}
