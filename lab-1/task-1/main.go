package task_1

import (
	"fmt"
	"math"
	"numerical-methods/linalg"
	"numerical-methods/utils"
	"os"
	"strings"
)

type Vector = linalg.Vector
type Matrix = linalg.Matrix

// region helpers
func maxRow(m Matrix, k int) int {
	pivot := k

	for i := k + 1; i < len(m); i++ {
		if math.Abs(m[i][k]) > math.Abs(m[pivot][k]) {
			pivot = i
		}
	}

	return pivot
}

//endregion helpers

func LUDecompose(matrix Matrix) (Matrix, Matrix, Matrix, int, error) {
	n := len(matrix)
	if n == 0 {
		return nil, nil, nil, 0, fmt.Errorf("matrix is empty")
	}

	a := matrix.Copy()
	p := linalg.IdentityMatrix(n)
	l := linalg.IdentityMatrix(n)
	swaps := 0

	for k := 0; k < n-1; k++ {
		pivot := maxRow(a, k)

		if math.Abs(a[pivot][k]) < utils.Eps {
			return nil, nil, nil, 0, fmt.Errorf("pivot is so small %d", k+1)
		}

		pStep := linalg.IdentityMatrix(n)
		pStep[k], pStep[pivot] = pStep[pivot], pStep[k]

		aHat := pStep.Multiply(a)

		p = pStep.Multiply(p)
		if pivot != k {
			swaps++
		}

		m := linalg.IdentityMatrix(n)
		mInverse := linalg.IdentityMatrix(n)

		for i := k + 1; i < n; i++ {
			mu := aHat[i][k] / aHat[k][k]

			m[i][k] = -mu
			mInverse[i][k] = mu
		}

		a = m.Multiply(aHat)

		l = pStep.Multiply(l).Multiply(pStep).Multiply(mInverse)
	}

	if math.Abs(a[n-1][n-1]) < utils.Eps {
		return nil, nil, nil, 0, fmt.Errorf("matrix is degenerate")
	}

	u := a.Copy()
	return p, l, u, swaps, nil
}

func SolveLowerTriangular(l Matrix, b Vector) (Vector, error) {
	z := make(Vector, len(l))

	for i := range l {
		if math.Abs(l[i][i]) < utils.Eps {
			return nil, fmt.Errorf("l[%d][%d] = 0", i+1, i+1)
		}

		value := b[i]
		for j := 0; j < i; j++ {
			value -= l[i][j] * z[j]
		}
		z[i] = value / l[i][i]
	}

	return z, nil
}

func SolveUpperTriangular(u Matrix, b Vector) (Vector, error) {
	x := make(Vector, len(u))

	for i := len(u) - 1; i >= 0; i-- {
		if math.Abs(u[i][i]) < utils.Eps {
			return nil, fmt.Errorf("u[%d][%d] = 0", i+1, i+1)
		}

		value := b[i]
		for j := i + 1; j < len(u); j++ {
			value -= u[i][j] * x[j]
		}
		x[i] = value / u[i][i]
	}

	return x, nil
}

func inverseMatrix(p, l, u Matrix) (Matrix, error) {
	inverse := linalg.ZeroMatrix(len(p))

	for column := 0; column < len(p); column++ {
		e := make(Vector, len(p))
		e[column] = 1

		z, _ := SolveLowerTriangular(l, p.MultiplyVector(e))
		x, err := SolveUpperTriangular(u, z)
		if err != nil {
			return nil, err
		}

		for row := range inverse {
			inverse[row][column] = x[row]
		}
	}

	return inverse, nil
}

func detLU(m Matrix, swaps int) float64 {
	determinant := 1.0 - 2.0*float64(swaps%2)

	for i := range m {
		determinant *= m[i][i]
	}

	return determinant
}

func maxIdentityError(a Matrix) float64 {
	maximum := 0.0

	for i := range a {
		for j := range a[i] {
			expected := 0.0
			if i == j {
				expected = 1
			}
			maximum = math.Max(maximum, math.Abs(a[i][j]-expected))
		}
	}

	return maximum
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
	a := linalg.ZeroMatrix(n)
	b := make(Vector, n)

	for i := 0; i < n; i++ {
		row, rowErr := utils.ParseRow(lines[i], n+1)
		if rowErr != nil {
			return fmt.Errorf("line %d of input file: %w", i+1, rowErr)
		}

		copy(a[i], row[:n])
		b[i] = row[n]
	}

	p, l, u, swaps, err := LUDecompose(a)
	if err != nil {
		return err
	}

	z, err := SolveLowerTriangular(l, p.MultiplyVector(b))
	if err != nil {
		return err
	}
	x, err := SolveUpperTriangular(u, z)
	if err != nil {
		return err
	}

	luProduct := l.Multiply(u)
	paProduct := p.Multiply(a)
	inverse, err := inverseMatrix(p, l, u)
	if err != nil {
		return err
	}
	aInverseProduct := a.Multiply(inverse)

	fmt.Println("Matrix P:")
	p.Print()
	fmt.Println("\nMatrix L:")
	l.Print()
	fmt.Println("\nMatrix U:")
	u.Print()
	fmt.Println("\nResult of LU:")
	luProduct.Print()
	fmt.Println("\nSolution of Ax = b:")
	x.Print()
	fmt.Println("\nDeterminant A:", fmt.Sprintf("%.8f", detLU(u, swaps)))
	fmt.Println("\nInverse matrix A^-1:")
	inverse.Print()
	fmt.Println("\nCheck of AA^-1 = E:")
	aInverseProduct.Print()
	inverseError := maxIdentityError(aInverseProduct)
	fmt.Println("max |AA^-1 − E| =", fmt.Sprintf("%.8f", inverseError))

	fmt.Println("\nCheck LU = PA:")
	fmt.Println("PA:")
	paProduct.Print()
	luError := luProduct.MaxDifference(paProduct)
	fmt.Println("max | LU − PA | =", fmt.Sprintf("%.8f", luError))
	fmt.Println("Check passed:", luError < utils.Eps)

	return nil
}
