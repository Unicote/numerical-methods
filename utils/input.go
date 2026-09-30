package utils

import (
	"fmt"
	"numerical-methods/linalg"
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
