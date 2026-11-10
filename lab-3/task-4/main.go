package task_4

import (
	"fmt"
	"math"
	"numerical-methods/utils"
)

func NearestNode(xs []float64, x float64) int {
	best := 0
	for i := range xs {
		if math.Abs(xs[i]-x) < math.Abs(xs[best]-x) {
			best = i
		}
	}
	return best
}

func Derivatives(xs, ys []float64, j int, x float64) (left, right, first, second float64) {
	left = (ys[j] - ys[j-1]) / (xs[j] - xs[j-1])
	right = (ys[j+1] - ys[j]) / (xs[j+1] - xs[j])
	change := (right - left) / (xs[j+1] - xs[j-1])
	first = left + change*(2*x-xs[j-1]-xs[j])
	second = 2 * change
	return
}

func Run(inputPath string) error {
	data, err := utils.ReadRows(inputPath, 5, 5, 1)
	if err != nil {
		return err
	}
	xs, ys := data[0], data[1]
	xStar := data[2][0]
	for i := 1; i < len(xs); i++ {
		if xs[i] <= xs[i-1] {
			return fmt.Errorf("nodes must be strictly increasing")
		}
	}
	if xStar < xs[0] || xStar > xs[len(xs)-1] {
		return fmt.Errorf("x* = %g is outside the table", xStar)
	}

	fmt.Println("Исходная таблица:")
	fmt.Printf("%3s  %8s  %10s\n", "i", "x_i", "y_i")
	for i := range xs {
		fmt.Printf("%3d  %8.3f  %10.5f\n", i, xs[i], ys[i])
	}
	fmt.Println()
	fmt.Printf("Точка вычисления: X* = %.4f\n", xStar)

	j := min(max(NearestNode(xs, xStar), 1), len(xs)-2)
	fmt.Printf("Использованные узлы: x_%d = %.3f, x_%d = %.3f, x_%d = %.3f\n", j-1, xs[j-1], j, xs[j], j+1, xs[j+1])
	fmt.Println()

	left, right, first, second := Derivatives(xs, ys, j, xStar)
	fmt.Println("Первая производная:")
	fmt.Printf("		левосторонняя (1-й порядок):   y'(X*) = (y_%d - y_%d)/(x_%d - x_%d) = %.8f\n", j, j-1, j, j-1, left)
	fmt.Printf("		правосторонняя (1-й порядок):  y'(X*) = (y_%d - y_%d)/(x_%d - x_%d) = %.8f\n", j+1, j, j+1, j, right)
	fmt.Printf("		2-й порядок точности:          y'(X*) = %.8f\n", first)
	fmt.Println()
	fmt.Println("Вторая производная:")
	fmt.Printf("		y''(X*) = 2*[(y_%d - y_%d)/(x_%d - x_%d) - (y_%d - y_%d)/(x_%d - x_%d)] / (x_%d - x_%d) = %.8f\n", j+1, j, j+1, j, j, j-1, j, j-1, j+1, j-1, second)
	return nil
}
