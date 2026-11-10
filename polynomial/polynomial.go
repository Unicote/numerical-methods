package polynomial

import (
	"fmt"
	"math"
	"numerical-methods/utils"
	"slices"
	"strings"
)

type Polynomial []float64

func (p Polynomial) Value(x float64) float64 {
	result := 0.0
	for _, pCoef := range slices.Backward(p) {
		result = result*x + pCoef
	}
	return result
}

func (p Polynomial) String() string {
	var b strings.Builder
	isFirst := true
	for i, c := range slices.Backward(p) {
		if math.Abs(c) < utils.Eps {
			continue
		}
		b.WriteString(SignedTerm(c, isFirst))
		switch {
		case i == 1:
			b.WriteString("*x")
		case i > 1:
			b.WriteString(fmt.Sprintf("*x^%d", i))
		}
		isFirst = false
	}
	if isFirst {
		return "0"
	}
	return b.String()
}

// FormatFactor returns the "x - a" string
func FormatFactor(a float64) string {
	switch {
	case math.Abs(a) < utils.Eps:
		return "x"
	case a < 0:
		return fmt.Sprintf("(x + %.2f)", -a)
	default:
		return fmt.Sprintf("(x - %.2f)", a)
	}
}

func SignedTerm(c float64, first bool) string {
	switch {
	case first && c < 0:
		return fmt.Sprintf("-%.6f", -c)
	case first:
		return fmt.Sprintf("%.6f", c)
	case c < 0:
		return fmt.Sprintf(" - %.6f", -c)
	default:
		return fmt.Sprintf(" + %.6f", c)
	}
}
