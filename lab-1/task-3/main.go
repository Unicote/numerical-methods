package task_3

import (
	"fmt"
	"math"
	"numerical-methods/lab-1/task-1"
	"numerical-methods/linalg"
	"numerical-methods/utils"
)

type Vector = linalg.Vector
type Matrix = linalg.Matrix

const maxIterations = math.MaxInt32

// region Helpers
func UpperPart(alpha Matrix) Matrix {
	c := linalg.ZeroMatrix(len(alpha))

	for i := range alpha {
		for j := i + 1; j < len(alpha); j++ {
			c[i][j] = alpha[i][j]
		}
	}

	return c
}

// endregion Helpers

func JacobiMethod(a Matrix, b Vector) (Matrix, Vector, error) {
	n := len(a)
	alpha := linalg.ZeroMatrix(n)
	beta := make(Vector, n)

	for i := 0; i < n; i++ {
		if math.Abs(a[i][i]) < utils.Eps {
			return nil, nil, fmt.Errorf("a[%d][%d] = 0", i+1, i+1)
		}

		beta[i] = b[i] / a[i][i]

		for j := 0; j < n; j++ {
			if i != j {
				alpha[i][j] = -a[i][j] / a[i][i]
			}
		}
	}

	return alpha, beta, nil
}

func JacobiError(alphaNorm, difference float64) float64 {
	return alphaNorm / (1 - alphaNorm) * difference
}

func SeidelError(alphaNorm, cNorm, difference float64) float64 {
	return cNorm / (1 - alphaNorm) * difference
}

func SimpleIterations(alpha Matrix, beta Vector, epsilon, alphaNorm float64) (Vector, int, float64, error) {
	n := len(alpha)
	previous := append(Vector(nil), beta...)

	for iteration := 1; iteration <= maxIterations; iteration++ {
		current := make(Vector, n)

		for i := 0; i < n; i++ {
			value := beta[i]

			for j := 0; j < n; j++ {
				value += alpha[i][j] * previous[j]
			}

			current[i] = value
		}

		difference := current.MaxDifference(previous)
		criterion := difference

		if alphaNorm < 1 {
			criterion = JacobiError(alphaNorm, difference)
		}

		if criterion <= epsilon {
			return current, iteration, criterion, nil
		}

		previous = current
	}

	return nil, 0, 0, fmt.Errorf("Simple iterations does not complete on %d ", maxIterations) // TODO
}

func Seidel(alpha Matrix, beta Vector, epsilon, alphaNorm, cNorm float64) (Vector, int, float64, error) {
	n := len(alpha)
	previous := append(Vector(nil), beta...)

	for iteration := 1; iteration <= maxIterations; iteration++ {
		current := make(Vector, n)

		for i := 0; i < n; i++ {
			value := beta[i]

			for j := 0; j < i; j++ {
				value += alpha[i][j] * current[j]
			}

			for j := i + 1; j < n; j++ {
				value += alpha[i][j] * previous[j]
			}

			current[i] = value
		}

		difference := current.MaxDifference(previous)
		criterion := difference

		if alphaNorm < 1 {
			criterion = SeidelError(alphaNorm, cNorm, difference)
		}

		if criterion <= epsilon {
			return current, iteration, criterion, nil
		}

		previous = current
	}

	return nil, 0, 0, fmt.Errorf("m Seidel does not complete on %d", maxIterations)
}

func LU(matrix Matrix, b Vector) (Vector, error) {
	p, l, u, _, err := task_1.LUDecompose(matrix)
	if err != nil {
		return nil, err
	}

	pb := p.MultiplyVector(b)
	z, err := task_1.SolveLowerTriangular(l, pb)
	if err != nil {
		return nil, err
	}

	return task_1.SolveUpperTriangular(u, z)
}

func Run(inputPath string) error {
	lines, err := utils.ReadLines(inputPath)
	if err != nil {
		return err
	}

	if len(lines) == 0 {
		return fmt.Errorf("input file is empty")
	}

	n := len(lines)
	a := linalg.ZeroMatrix(n)
	b := make(Vector, n)

	for i, line := range lines {
		row, rowErr := utils.ParseRow(line, n+1)
		if rowErr != nil {
			return fmt.Errorf("line %d of input file: %w", i+1, rowErr)
		}

		copy(a[i], row[:n])
		b[i] = row[n]
	}

	alpha, beta, err := JacobiMethod(a, b)
	if err != nil {
		return err
	}
	alphaNorm := alpha.Norm()
	c := UpperPart(alpha)
	cNorm := c.Norm()

	refSolution, err := LU(a, b)
	if err != nil {
		return err
	}

	simple, simpleIterations, simpleError, err := SimpleIterations(alpha, beta, utils.Eps, alphaNorm)
	if err != nil {
		return err
	}

	seidel, seidelIterations, seidelError, err := Seidel(alpha, beta, utils.Eps, alphaNorm, cNorm)
	if err != nil {
		return err
	}

	fmt.Println("Epsilon:", utils.Eps)
	fmt.Println("\nMatrix alpha:")
	alpha.Print()
	fmt.Println("\nVector beta:")
	beta.Print()
	fmt.Println("Norm ||alpha||:", alphaNorm)

	fmt.Println("\nReference solution:")
	refSolution.Print()

	fmt.Println("\nSimple iterations method:")
	simple.Print()
	fmt.Println("Iterations:", simpleIterations)
	if alphaNorm < 1 {
		fmt.Println("Error estimate:", simpleError)
	} else {
		fmt.Println("max |x^(k) - x^(k-1)|:", simpleError)
	}
	fmt.Println("max |x - xReference|:", simple.MaxDifference(refSolution))

	fmt.Println("\nSeidel method:")
	seidel.Print()
	fmt.Println("Iterations:", seidelIterations)
	if alphaNorm < 1 {
		fmt.Println("Error estimate:", seidelError)
	} else {
		fmt.Println("max |x^(k) - x^(k-1)|:", seidelError)
	}
	fmt.Println("max |x - xReference|:", seidel.MaxDifference(refSolution))

	fmt.Println("\nConvergence speed comparison:")
	if simpleIterations < seidelIterations {
		fmt.Println("Simple iterations method converged faster")
	} else if seidelIterations < simpleIterations {
		fmt.Println("Seidel method converged faster")
	} else {
		fmt.Println("Both methods used the same number of iterations")
	}

	return nil
}
