package task_4

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
func IsSymmetric(a Matrix) bool {
	for i := range a {
		for j := i + 1; j < len(a); j++ {
			if math.Abs(a[i][j]-a[j][i]) > utils.Eps {
				return false
			}
		}
	}

	return true
}

// endregion helpers

func t(a Matrix) float64 {
	sum := 0.0

	for i := range a {
		for j := i + 1; j < len(a); j++ {
			sum += a[i][j] * a[i][j]
		}
	}

	return math.Sqrt(sum)
}

func MaxOffDiagonalIndex(a Matrix) (int, int) {
	row := 0
	column := 1
	maximum := math.Abs(a[row][column])

	for i := range a {
		for j := i + 1; j < len(a); j++ {
			if math.Abs(a[i][j]) > maximum {
				maximum = math.Abs(a[i][j])
				row = i
				column = j
			}
		}
	}

	return row, column
}

func Phi(a Matrix, i, j int) float64 {
	diff := a[i][i] - a[j][j]

	if math.Abs(diff) < utils.Eps {
		if a[i][j] < 0 {
			return -math.Pi / 4
		}
		return math.Pi / 4
	}

	return 0.5 * math.Atan(2*a[i][j]/diff)
}

func U(n, i, j int, phi float64) Matrix {
	u := linalg.IdentityMatrix(n)
	cos := math.Cos(phi)
	sin := math.Sin(phi)

	u[i][i] = cos
	u[i][j] = -sin
	u[j][i] = sin
	u[j][j] = cos

	return u
}

func DiagonalMatrix(values Vector) Matrix {
	lambda := linalg.ZeroMatrix(len(values))

	for i := range values {
		lambda[i][i] = values[i]
	}

	return lambda
}

func JacobiRotations(a Matrix, epsilon float64) (Vector, Matrix, int, error) {
	if len(a) == 0 {
		return nil, nil, 0, fmt.Errorf("matrix is empty")
	}
	if !IsSymmetric(a) {
		return nil, nil, 0, fmt.Errorf("Jacobi rotations method applies only to symmetric matrices")
	}

	current := a.Copy()
	eigenvectors := linalg.IdentityMatrix(len(a))

	for iteration := 0; iteration < maxIterations; iteration++ {
		if t(current) <= epsilon {
			eigenvalues := make(Vector, len(current))
			for i := range current {
				eigenvalues[i] = current[i][i]
			}

			return eigenvalues, eigenvectors, iteration, nil
		}

		i, j := MaxOffDiagonalIndex(current)
		phi := Phi(current, i, j)
		u := U(len(a), i, j, phi)

		current = u.Transpose().Multiply(current).Multiply(u)
		eigenvectors = eigenvectors.Multiply(u)
	}

	return nil, nil, 0, fmt.Errorf("Jacobi rotations did not finish after %d iterations", maxIterations)
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

	eigenvalues, eigenvectors, iterations, err := JacobiRotations(a, utils.Eps)
	if err != nil {
		return err
	}

	lambda := DiagonalMatrix(eigenvalues)
	leftSide := a.Multiply(eigenvectors)
	rightSide := eigenvectors.Multiply(lambda)
	checkError := leftSide.MaxDifference(rightSide)

	fmt.Println("Epsilon:", utils.Eps)
	fmt.Println("Iterations:", iterations)

	fmt.Println("\nEigenvalues:")
	eigenvalues.Print()
	fmt.Println("\nDiagonal matrix Lambda:")
	lambda.Print()

	fmt.Println("\nEigenvectors:")
	for column := 0; column < n; column++ {
		vector := make(Vector, n)
		for row := 0; row < n; row++ {
			vector[row] = eigenvectors[row][column]
		}

		fmt.Printf("x%d = ", column+1)
		vector.Print()
	}

	fmt.Println("\nEigenvectors matrix V:")
	eigenvectors.Print()

	fmt.Println("\nCheck AV = V Lambda:")
	fmt.Println("AV:")
	leftSide.Print()
	fmt.Println("V * Lambda:")
	rightSide.Print()
	fmt.Println("max |A V - V * Lambda|:", checkError)

	return nil
}
