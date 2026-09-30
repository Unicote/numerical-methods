package linalg

import "math"

type Matrix []Vector

func ZeroMatrix(n int) Matrix {
	matrix := make(Matrix, n)
	for i := range matrix {
		matrix[i] = make(Vector, n)
	}
	return matrix
}

func IdentityMatrix(n int) Matrix {
	matrix := ZeroMatrix(n)
	for i := 0; i < n; i++ {
		matrix[i][i] = 1
	}
	return matrix
}

func (m Matrix) Copy() Matrix {
	copyMatrix := make(Matrix, len(m))
	for i := range m {
		copyMatrix[i] = append(Vector(nil), m[i]...)
	}
	return copyMatrix
}

func (m Matrix) Print() {
	for _, row := range m {
		row.Print()
	}
}

func (m Matrix) Multiply(b Matrix) Matrix {
	n := len(m)
	result := ZeroMatrix(n)

	for i := 0; i < n; i++ {
		for k := 0; k < n; k++ {
			for j := 0; j < n; j++ {
				result[i][j] += m[i][k] * b[k][j]
			}
		}
	}

	return result
}

func (m Matrix) MultiplyVector(v Vector) Vector {
	result := make(Vector, len(m))

	for i := range m {
		for j := range m {
			result[i] += m[i][j] * v[j]
		}
	}

	return result
}

func (m Matrix) Norm() float64 {
	maximum := 0.0

	for i := range m {
		sum := 0.0
		for j := range m[i] {
			sum += math.Abs(m[i][j])
		}
		maximum = math.Max(maximum, sum)
	}

	return maximum
}

func (m Matrix) MaxDifference(second Matrix) float64 {
	maximum := 0.0

	for i := range m {
		for j := range m[i] {
			maximum = math.Max(maximum, math.Abs(m[i][j]-second[i][j]))
		}
	}

	return maximum
}

func (m Matrix) Transpose() Matrix {
	transposed := ZeroMatrix(len(m))

	for i := range m {
		for j := range m[i] {
			transposed[j][i] = m[i][j]
		}
	}

	return transposed
}
