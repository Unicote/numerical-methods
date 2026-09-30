package task_5

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

const maxIterations = math.MaxInt32

// region helpers

func IsSquare(a Matrix) bool {
	if len(a) == 0 {
		return false
	}

	for i := range a {
		if len(a[i]) != len(a) {
			return false
		}
	}

	return true
}

func ComplexMaxDifference(first, second []complex128) float64 {
	maximum := 0.0

	for i := range first {
		maximum = math.Max(maximum, complexAbs(first[i]-second[i]))
	}

	return maximum
}

func complexAbs(value complex128) float64 {
	return math.Hypot(real(value), imag(value))
}

func PrintEigenvalues(eigenvalues []complex128) {
	fmt.Print("[")
	for i, value := range eigenvalues {
		if i > 0 {
			fmt.Print(", ")
		}

		if math.Abs(imag(value)) < utils.Eps {
			fmt.Printf("%.8f", real(value))
		} else {
			fmt.Printf("%.8f%+.8fi", real(value), imag(value))
		}
	}
	fmt.Println("]")
}

// endregion helpers

func HouseholderMatrix(a Matrix, column int) Matrix {
	n := len(a)
	v := make(Vector, n)
	bNorm := 0.0

	for i := column; i < n; i++ {
		bNorm += a[i][column] * a[i][column]
	}
	bNorm = math.Sqrt(bNorm)

	if bNorm < utils.Eps {
		return linalg.IdentityMatrix(n)
	}

	sign := 1.0
	if a[column][column] < 0 {
		sign = -1
	}

	v[column] = a[column][column] + sign*bNorm
	for i := column + 1; i < n; i++ {
		v[i] = a[i][column]
	}

	denominator := 0.0
	for i := column; i < n; i++ {
		denominator += v[i] * v[i]
	}

	if denominator < utils.Eps {
		return linalg.IdentityMatrix(n)
	}

	h := linalg.IdentityMatrix(n)
	for i := column; i < n; i++ {
		for j := column; j < n; j++ {
			h[i][j] -= 2 * v[i] * v[j] / denominator
		}
	}

	return h
}

func QRDecompose(a Matrix) (Matrix, Matrix, error) {
	if !IsSquare(a) {
		return nil, nil, fmt.Errorf("matrix must be square")
	}

	r := a.Copy()
	q := linalg.IdentityMatrix(len(a))

	for column := 0; column < len(a)-1; column++ {
		h := HouseholderMatrix(r, column)
		r = h.Multiply(r)
		q = q.Multiply(h)
	}

	for i := range r {
		for j := 0; j < i; j++ {
			if math.Abs(r[i][j]) < utils.Eps {
				r[i][j] = 0
			}
		}
	}

	return q, r, nil
}

func Eigenvalues2x2(a, b, c, d float64) (complex128, complex128) {
	trace := a + d
	determinant := a*d - b*c
	discriminant := trace*trace - 4*determinant

	if discriminant >= 0 {
		squareRoot := math.Sqrt(discriminant)
		return complex((trace+squareRoot)/2, 0), complex((trace-squareRoot)/2, 0)
	}

	realPart := trace / 2
	imaginaryPart := math.Sqrt(-discriminant) / 2
	return complex(realPart, imaginaryPart), complex(realPart, -imaginaryPart)
}

func GetEigenvalues(a Matrix, epsilon float64) ([]complex128, bool, bool) {
	eigenvalues := make([]complex128, 0, len(a))
	hasComplexPair := false

	for i := 0; i < len(a); {
		if i == len(a)-1 || math.Abs(a[i+1][i]) <= epsilon {
			eigenvalues = append(eigenvalues, complex(a[i][i], 0))
			i++
			continue
		}

		if i+2 < len(a) && math.Abs(a[i+2][i+1]) > epsilon {
			return nil, false, false
		}

		first, second := Eigenvalues2x2(a[i][i], a[i][i+1], a[i+1][i], a[i+1][i+1])
		if math.Abs(imag(first)) <= epsilon && math.Abs(imag(second)) <= epsilon {
			return nil, false, false
		}

		eigenvalues = append(eigenvalues, first, second)
		hasComplexPair = true
		i += 2
	}

	return eigenvalues, true, hasComplexPair
}

func QREigenvalues(a Matrix, epsilon float64) ([]complex128, int, error) {
	if !IsSquare(a) {
		return nil, 0, fmt.Errorf("matrix must be square")
	}

	current := a.Copy()
	var previousComplexEigenvalues []complex128

	for iteration := 0; iteration < maxIterations; iteration++ {
		q, r, err := QRDecompose(current)
		if err != nil {
			return nil, 0, err
		}

		current = r.Multiply(q) // a^(k+1)= r^k q^k
		eigenvalues, completed, hasComplexPair := GetEigenvalues(current, epsilon)
		if !completed {
			continue
		}

		if !hasComplexPair {
			return eigenvalues, iteration + 1, nil
		}

		if len(previousComplexEigenvalues) == len(eigenvalues) && ComplexMaxDifference(eigenvalues, previousComplexEigenvalues) <= epsilon {
			return eigenvalues, iteration + 1, nil
		}

		previousComplexEigenvalues = eigenvalues
	}

	return nil, 0, fmt.Errorf("QR algorithm did not finish after %d iterations", maxIterations)
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

	for i, line := range lines {
		row, rowErr := utils.ParseRow(line, n)
		if rowErr != nil {
			return fmt.Errorf("line %d of input file: %w", i+1, rowErr)
		}

		copy(a[i], row)
	}

	q, r, err := QRDecompose(a)
	if err != nil {
		return err
	}

	eigenvalues, iterations, err := QREigenvalues(a, utils.Eps)
	if err != nil {
		return err
	}

	fmt.Println("Input matrix A:")
	a.Print()

	fmt.Println("\nMatrix Q:")
	q.Print()

	fmt.Println("\nMatrix R:")
	r.Print()

	fmt.Println("\nResult of QR:")
	q.Multiply(r).Print()
	fmt.Println("max |A - QR|:", a.MaxDifference(q.Multiply(r)))

	fmt.Println("\nEigenvalues:")
	PrintEigenvalues(eigenvalues)
	fmt.Println("QR iterations:", iterations)

	return nil
}
