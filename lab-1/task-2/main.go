package task_2

import (
	"fmt"
	"math"
	"numerical-methods/linalg"
	"numerical-methods/utils"
	"os"
	"strings"
)

type Vector = linalg.Vector

func Solve(a, b, c, d Vector) (Vector, Vector, Vector, error) {
	n := len(b)

	if n == 0 || len(a) != n || len(c) != n || len(d) != n {
		return nil, nil, nil, fmt.Errorf("incorrect number of diagonal values")
	}

	if math.Abs(a[0]) >= utils.Eps || math.Abs(c[n-1]) >= utils.Eps {
		return nil, nil, nil, fmt.Errorf("values must be: lower[0] = 0; upper[n-1] = 0")
	}

	p := make(Vector, n)
	q := make(Vector, n)

	if math.Abs(b[0]) < utils.Eps {
		return nil, nil, nil, fmt.Errorf("division by null")
	}
	p[0] = -c[0] / b[0]
	q[0] = d[0] / b[0]

	for i := 1; i < n; i++ {
		denominator := a[i]*p[i-1] + b[i]

		if math.Abs(denominator) < utils.Eps {
			return nil, nil, nil, fmt.Errorf("division by null on %d", i+1)
		}

		p[i] = -c[i] / denominator
		q[i] = (d[i] - a[i]*q[i-1]) / denominator
	}

	x := make(Vector, n)
	x[n-1] = q[n-1]

	for i := n - 2; i >= 0; i-- {
		x[i] = p[i]*x[i+1] + q[i]
	}

	return p, q, x, nil
}

func Run(inputPath string) error {
	data, err := os.ReadFile(inputPath)
	if err != nil {
		return fmt.Errorf("read input file: %w", err)
	}

	lines := make([]string, 0)
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimSpace(line)
		if line != "" {
			lines = append(lines, line)
		}
	}

	if len(lines) == 0 {
		return fmt.Errorf("input file is empty")
	}

	n := len(lines)
	a := make(Vector, n) // lower diagonal
	b := make(Vector, n) //middle diagonal
	c := make(Vector, n) // upper diagonal
	d := make(Vector, n)

	const countVar = 4

	for i, line := range lines {
		row, rowErr := utils.ParseRow(line, countVar)
		if rowErr != nil {
			return fmt.Errorf("line %d of input file: %w", i+1, rowErr)
		}

		a[i] = row[0]
		b[i] = row[1]
		c[i] = row[2]
		d[i] = row[3]
	}

	p, q, x, err := Solve(a, b, c, d)
	if err != nil {
		return err
	}

	fmt.Println("P coefficient:")
	p.Print()

	fmt.Println("Q coefficient:")
	q.Print()

	fmt.Println("Solution of system:")
	x.Print()

	return nil
}
