package linalg

import (
	"fmt"
	"math"
)

type Vector []float64

func (v Vector) Print() {
	fmt.Print("[")
	for i, value := range v {
		if i > 0 {
			fmt.Print(", ")
		}
		fmt.Print(fmt.Sprintf("%.8f", value))
	}
	fmt.Println("]")
}

func (v Vector) MaxDifference(second Vector) float64 {
	maximum := 0.0

	for i := range v {
		maximum = math.Max(maximum, math.Abs(v[i]-second[i]))
	}

	return maximum
}
