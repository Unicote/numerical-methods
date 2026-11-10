package task_3

import (
	"fmt"
	"math"
	lu "numerical-methods/lab-1/task-3"
	"numerical-methods/linalg"
	"numerical-methods/plot"
	"numerical-methods/polynomial"
	"numerical-methods/utils"
	"path/filepath"
)

// region helpers
func printSystem(matrix linalg.Matrix, right linalg.Vector) {
	for i := range matrix {
		fmt.Print("  ")
		for j := range matrix[i] {
			fmt.Printf("%12.5f*a%d", matrix[i][j], j)
			if j < len(matrix[i])-1 {
				fmt.Print(" +")
			}
		}
		fmt.Printf(" = %12.6f\n", right[i])
	}
}

// endregion helpers

type Polynomial = polynomial.Polynomial

func System(xs, ys []float64, n int) (linalg.Matrix, linalg.Vector) {
	size := n + 1
	matrix := linalg.ZeroMatrix(size)
	right := make(linalg.Vector, size)
	for i := range size {
		for j := range size {
			for _, x := range xs {
				matrix[i][j] += math.Pow(x, float64(i+j))
			}
		}
		for k, x := range xs {
			right[i] += ys[k] * math.Pow(x, float64(i))
		}
	}
	return matrix, right
}

func SolveSystem(xs, ys []float64, degree int) (Polynomial, linalg.Matrix, linalg.Vector, error) {
	matrix, right := System(xs, ys, degree)
	coefficients, err := lu.LU(matrix.Copy(), right)
	if err != nil {
		return nil, nil, nil, err
	}
	return Polynomial(coefficients), matrix, right, nil
}

func SquaredErrorSum(p Polynomial, xs, ys []float64) float64 {
	sum := 0.0
	for k := range xs {
		diff := p.Value(xs[k]) - ys[k]
		sum += diff * diff
	}
	return sum
}

func Run(inputPath string) error {
	data, err := utils.ReadRows(inputPath, 6, 6)
	if err != nil {
		return err
	}
	xs, ys := data[0], data[1]

	fmt.Println("Исходные данные:")
	fmt.Printf("%3s  %10s  %12s\n", "i", "x_i", "y_i")
	for i := range xs {
		fmt.Printf("%3d  %10.4f  %12.6f\n", i, xs[i], ys[i])
	}
	fmt.Println()

	degrees := []int{1, 2, 3}
	polynomials := make([]Polynomial, len(degrees))
	squaredSums := make([]float64, len(degrees))
	mse := make([]float64, len(degrees))

	for index, degree := range degrees {
		p, matrix, right, err := SolveSystem(xs, ys, degree)
		if err != nil {
			return fmt.Errorf("degree %d: %w", degree, err)
		}
		polynomials[index] = p
		squaredSums[index] = SquaredErrorSum(p, xs, ys)
		mse[index] = squaredSums[index] / float64(len(xs))

		fmt.Printf("---- Многочлен степени %d ----\n", degree)
		fmt.Println("Нормальная система:")
		printSystem(matrix, right)
		fmt.Println("Коэффициенты:")
		for j, c := range p {
			fmt.Printf("  a%d = %.8f\n", j, c)
		}
		fmt.Printf("P%d(x) = %s\n", degree, p.String())
		fmt.Printf("Сумма квадратов ошибок Ф%d = %.8f\n", degree, squaredSums[index])
		fmt.Printf("MSE%d = Ф%d / N = %.8f\n", degree, degree, mse[index])
		fmt.Println()
	}

	fmt.Println("Сравнение ошибок:")
	fmt.Printf("%9s  %14s  %14s\n", "Степень", "Ф", "MSE")
	best := 0
	for index, degree := range degrees {
		fmt.Printf("%9d  %14.8f  %14.8f\n", degree, squaredSums[index], mse[index])
		if mse[index] < mse[best] {
			best = index
		}
	}
	fmt.Printf("Лучше всего данные приближает P%d(x): MSE%d = %.8f\n", degrees[best], degrees[best], mse[best])
	fmt.Println()

	a, b := xs[0], xs[len(xs)-1]
	margin := 0.05 * (b - a)
	figure := plot.Figure{Title: "Метод наименьших квадратов", XLabel: "x", YLabel: "y"}
	colors := []string{"#1f77b4", "#2ca02c", "#d62728"}
	for index, degree := range degrees {
		px, py := plot.Sample(polynomials[index].Value, a-margin, b+margin, 300)
		figure.Series = append(figure.Series, plot.Series{
			Name: fmt.Sprintf("P%d(x), MSE=%.5f", degree, mse[index]), X: px, Y: py, Color: colors[index], Width: 2.5,
		})
	}
	figure.Series = append(figure.Series, plot.Series{Name: "Исходные точки", X: xs, Y: ys, Points: true, Color: "#222222", Width: 6})
	plotPath := filepath.Join(filepath.Dir(inputPath), "plot.svg")
	if err := figure.Save(plotPath); err != nil {
		return fmt.Errorf("save plot: %w", err)
	}
	return nil
}
