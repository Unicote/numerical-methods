package utils

import (
	"fmt"
	"math"
	"strings"
	"unicode/utf8"
)

func PadRight(s string, width int) string {
	if n := utf8.RuneCountInString(s); n < width {
		return s + strings.Repeat(" ", width-n)
	}
	return s
}

func PadLeft(s string, width int) string {
	if n := utf8.RuneCountInString(s); n < width {
		return strings.Repeat(" ", width-n) + s
	}
	return s
}

func FormatOptional(value float64) string {
	if math.IsNaN(value) {
		return "-"
	}
	return fmt.Sprintf("%.6f", value)
}
