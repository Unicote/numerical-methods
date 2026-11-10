package task_1

import (
	"fmt"
	"math"
	"numerical-methods/plot"
	"numerical-methods/polynomial"
	"numerical-methods/utils"
	"path/filepath"
	"slices"
	"strings"
)

func Function(x float64) float64 {
	return math.Sin(x)
}

// region helpers
func printDividedDifferences(xs []float64, table [][]float64) {
	n := len(xs)
	fmt.Printf("%3s  %12s  %12s", "i", "x_i", "f(x_i)")
	for k := 1; k < n; k++ {
		fmt.Printf("  %s", utils.PadLeft(fmt.Sprintf("%d-го порядка", k), 14))
	}
	fmt.Println()
	for i := range n {
		fmt.Printf("%3d  %12.8f  %12.8f", i, xs[i], table[0][i])
		for k := 1; k < n-i; k++ {
			fmt.Printf("  %14.8f", table[k][i])
		}
		fmt.Println()
	}
}

//endregion helpers

// region Lagrange

func LagrangeWeights(xs, ys []float64) ([]float64, error) {
	result := make([]float64, len(xs))
	for i := range xs {
		denominator := 1.0
		for j := range xs {
			if j != i {
				denominator *= xs[i] - xs[j]
			}
		}
		if math.Abs(denominator) < utils.Eps {
			return nil, fmt.Errorf("interpolation nodes must be different")
		}
		result[i] = ys[i] / denominator
	}
	return result, nil
}

func LagrangeValue(xs, weights []float64, x float64) float64 {
	result := 0.0
	for i := range xs {
		term := weights[i]
		for j := range xs {
			if j != i {
				term *= x - xs[j]
			}
		}
		result += term
	}
	return result
}

func LagrangeString(xs, weights []float64) string {
	var b strings.Builder
	for i := range xs {
		b.WriteString(polynomial.SignedTerm(weights[i], i == 0))
		for j := range xs {
			if j != i {
				b.WriteString("*" + polynomial.FormatFactor(xs[j]))
			}
		}
	}
	return b.String()
}

//endregion Lagrange

// region Newton

func DividedDifferences(xs, ys []float64) ([][]float64, error) {
	n := len(xs)
	table := make([][]float64, n)
	table[0] = append([]float64(nil), ys...)
	for k := 1; k < n; k++ {
		table[k] = make([]float64, n-k)
		for i := 0; i < n-k; i++ {
			denominator := xs[i+k] - xs[i]
			if math.Abs(denominator) < utils.Eps {
				return nil, fmt.Errorf("interpolation nodes must be different")
			}
			table[k][i] = (table[k-1][i+1] - table[k-1][i]) / denominator
		}
	}
	return table, nil
}

func NewtonValue(xs []float64, table [][]float64, x float64) float64 {
	result := 0.0
	for k, x0 := range slices.Backward(xs) {
		result = result*(x-x0) + table[k][0]
	}
	return result
}

func NewtonString(xs []float64, table [][]float64) string {
	var b strings.Builder
	for k := range xs {
		b.WriteString(polynomial.SignedTerm(table[k][0], k == 0))
		for j := range k {
			b.WriteString("*" + polynomial.FormatFactor(xs[j]))
		}
	}
	return b.String()
}

//endregion Newton

func solveSet(name string, xs []float64, xStar float64, plotPath string) error {
	ys := make([]float64, len(xs))
	for i := range xs {
		ys[i] = Function(xs[i])
	}

	fmt.Printf("%s)\n", name)
	fmt.Println("Узлы интерполяции:")
	fmt.Printf("%3s  %12s  %12s  %12s\n", "i", "x_i", "x_i / pi", "y_i")
	for i := range xs {
		fmt.Printf("%3d  %12.8f  %12.8f  %12.8f\n", i, xs[i], xs[i]/math.Pi, ys[i])
	}
	fmt.Println()

	weights, err := LagrangeWeights(xs, ys)
	if err != nil {
		return err
	}
	fmt.Println("Интерполяционный многочлен Лагранжа:")
	fmt.Println("L(x) = " + LagrangeString(xs, weights))
	fmt.Println()

	table, err := DividedDifferences(xs, ys)
	if err != nil {
		return err
	}
	fmt.Println("Таблица разделённых разностей:")
	printDividedDifferences(xs, table)
	fmt.Println()
	fmt.Println("Интерполяционный многочлен Ньютона:")
	fmt.Println("P(x) = " + NewtonString(xs, table))
	fmt.Println()

	lagrangeValue := LagrangeValue(xs, weights, xStar)
	newtonValue := NewtonValue(xs, table, xStar)
	exact := Function(xStar)
	fmt.Printf("X* = %.8f (= %.4f*pi)\n", xStar, xStar/math.Pi)
	fmt.Printf("L(X*)                  = %.10f\n", lagrangeValue)
	fmt.Printf("P(X*)                  = %.10f\n", newtonValue)
	fmt.Printf("f(X*) точное решение   = %.10f\n", exact)
	fmt.Printf("|L(X*) - f(X*)|        = %.10e\n", math.Abs(lagrangeValue-exact))
	fmt.Printf("|P(X*) - f(X*)|        = %.10e\n", math.Abs(newtonValue-exact))
	fmt.Printf("|L(X*) - P(X*)|        = %.10e\n", math.Abs(lagrangeValue-newtonValue))
	fmt.Println()

	a := math.Min(xs[0], xStar)
	b := math.Max(xs[len(xs)-1], xStar)
	margin := 0.5 * (b - a)
	a, b = a-margin, b+margin

	fx, fy := plot.Sample(Function, a, b, 400)
	lx, ly := plot.Sample(func(x float64) float64 { return LagrangeValue(xs, weights, x) }, a, b, 400)
	nx, ny := plot.Sample(func(x float64) float64 { return NewtonValue(xs, table, x) }, a, b, 400)
	figure := plot.Figure{
		Title:  fmt.Sprintf("Интерполяция y = sin(x), набор %s", name),
		XLabel: "x",
		YLabel: "y",
		Series: []plot.Series{
			{Name: "f(x) = sin(x)", X: fx, Y: fy, Color: "#1f77b4", Width: 5},
			{Name: "Лагранж L(x)", X: lx, Y: ly, Color: "#d62728", Width: 2.5},
			{Name: "Ньютон P(x)", X: nx, Y: ny, Color: "#2ca02c", Dashed: true, Width: 2.5},
			{Name: "Узлы (x_i, y_i)", X: xs, Y: ys, Points: true, Color: "#222222"},
			{Name: "X* (L(X*))", X: []float64{xStar}, Y: []float64{lagrangeValue}, Points: true, Color: "#ff7f0e", Width: 7},
		},
	}
	if err := figure.Save(plotPath); err != nil {
		return fmt.Errorf("save plot: %w", err)
	}
	return nil
}

func Run(inputPath string) error {
	data, err := utils.ReadRows(inputPath, 4, 4, 1)
	if err != nil {
		return err
	}

	fmt.Println("Function: y = sin(x)")
	fmt.Println()

	directory := filepath.Dir(inputPath)
	if err := solveSet("а", data[0], data[2][0], filepath.Join(directory, "plot-a.svg")); err != nil {
		return err
	}
	return solveSet("б", data[1], data[2][0], filepath.Join(directory, "plot-b.svg"))
}
