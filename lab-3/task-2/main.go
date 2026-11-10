package task_2

import (
	"fmt"
	"math"
	tridiagonal "numerical-methods/lab-1/task-2"
	"numerical-methods/linalg"
	"numerical-methods/plot"
	"numerical-methods/utils"
	"path/filepath"
)

type Spline struct {
	X          []float64
	A, B, C, D []float64
}

func makeSpline(xs, fs []float64) (Spline, error) {
	n := len(xs) - 1
	if n < 2 {
		return Spline{}, fmt.Errorf("at least three nodes are required")
	}
	h := make([]float64, n+1)
	for i := 1; i <= n; i++ {
		h[i] = xs[i] - xs[i-1]
		if h[i] <= 0 {
			return Spline{}, fmt.Errorf("nodes must be strictly increasing")
		}
	}

	m := n - 1
	lower := make(linalg.Vector, m)
	main := make(linalg.Vector, m)
	upper := make(linalg.Vector, m)
	right := make(linalg.Vector, m)
	for k := range m {
		i := k + 2
		if k > 0 {
			lower[k] = h[i-1]
		}
		main[k] = 2 * (h[i-1] + h[i])
		if k < m-1 {
			upper[k] = h[i]
		}
		right[k] = 3 * ((fs[i]-fs[i-1])/h[i] - (fs[i-1]-fs[i-2])/h[i-1])
	}
	_, _, solution, err := tridiagonal.Solve(lower, main, upper, right)
	if err != nil {
		return Spline{}, err
	}

	c := make([]float64, n+2)
	for k := range m {
		c[k+2] = solution[k]
	}

	s := Spline{
		X: xs,
		A: make([]float64, n+1),
		B: make([]float64, n+1),
		C: make([]float64, n+1),
		D: make([]float64, n+1),
	}
	for i := 1; i <= n; i++ {
		s.A[i] = fs[i-1]
		s.B[i] = (fs[i]-fs[i-1])/h[i] - h[i]*(c[i+1]+2*c[i])/3
		s.C[i] = c[i]
		s.D[i] = (c[i+1] - c[i]) / (3 * h[i])
	}
	return s, nil
}

func (s Spline) Segments() int { return len(s.X) - 1 }

func (s Spline) Segment(x float64) int {
	n := s.Segments()
	for i := 1; i <= n; i++ {
		if x <= s.X[i] {
			return i
		}
	}
	return n
}

func (s Spline) ValueOnSegment(i int, x float64) float64 {
	t := x - s.X[i-1]
	return s.A[i] + s.B[i]*t + s.C[i]*t*t + s.D[i]*t*t*t
}

func (s Spline) SecondDerivativeOn(i int, x float64) float64 {
	t := x - s.X[i-1]
	return 2*s.C[i] + 6*s.D[i]*t
}

func (s Spline) Value(x float64) float64 {
	return s.ValueOnSegment(s.Segment(x), x)
}

func Run(inputPath string) error {
	data, err := utils.ReadRows(inputPath, 5, 5, 1)
	if err != nil {
		return err
	}
	xs, fs := data[0], data[1]
	xStar := data[2][0]
	if xStar < xs[0] || xStar > xs[len(xs)-1] {
		return fmt.Errorf("point x* = %g is outside [%g, %g]", xStar, xs[0], xs[len(xs)-1])
	}

	spline, err := makeSpline(xs, fs)
	if err != nil {
		return err
	}
	n := spline.Segments()

	fmt.Println("Исходные узлы:")
	fmt.Printf("%3s  %10s  %12s\n", "i", "x_i", "f_i")
	for i := range xs {
		fmt.Printf("%3d  %10.4f  %12.6f\n", i, xs[i], fs[i])
	}
	fmt.Println()
	fmt.Println("Коэффициенты сплайна:")
	fmt.Printf("%3s  %18s  %12s  %12s  %12s  %12s\n", "i", "[x_{i-1}, x_i]", "a_i", "b_i", "c_i", "d_i")
	for i := 1; i <= n; i++ {
		fmt.Printf("%3d  %18s  %12.6f  %12.6f  %12.6f  %12.6f\n", i, fmt.Sprintf("[%.2f, %.2f]", xs[i-1], xs[i]), spline.A[i], spline.B[i], spline.C[i], spline.D[i])
	}
	fmt.Println()

	fmt.Println("Многочлены на участках:")
	for i := 1; i <= n; i++ {
		fmt.Printf("S_%d(x) = %.6f %+.6f(x - %.2f) %+.6f(x - %.2f)^2 %+.6f(x - %.2f)^3,  x in [%.2f, %.2f]\n", i, spline.A[i], spline.B[i], xs[i-1], spline.C[i], xs[i-1], spline.D[i], xs[i-1], xs[i-1], xs[i])
	}
	fmt.Println()

	segment := spline.Segment(xStar)
	value := spline.Value(xStar)
	fmt.Printf("X* = %.6f\n", xStar)
	fmt.Printf("X* принадлежит участку i = %d: [%.2f, %.2f]\n", segment, xs[segment-1], xs[segment])
	fmt.Printf("S(X*) = %.10f\n", value)
	fmt.Println()

	fmt.Println("Проверка прохождения через узлы:")
	fmt.Printf("%3s  %10s  %12s  %14s  %14s  %12s\n", "i", "x_i", "f_i", "S_left(x_i)", "S_right(x_i)", "max |S - f|")
	for i := range xs {
		left, right := math.NaN(), math.NaN()
		maxError := 0.0
		if i > 0 {
			left = spline.ValueOnSegment(i, xs[i])
			maxError = math.Max(maxError, math.Abs(left-fs[i]))
		}
		if i < n {
			right = spline.ValueOnSegment(i+1, xs[i])
			maxError = math.Max(maxError, math.Abs(right-fs[i]))
		}
		fmt.Printf("%3d  %10.4f  %12.6f  %14s  %14s  %12.2e\n", i, xs[i], fs[i], utils.FormatOptional(left), utils.FormatOptional(right), maxError)
	}
	fmt.Println()

	fmt.Println("Проверка граничных условий:")
	fmt.Printf("S''(x_0 = %.2f) = 2c_1 = %.10f\n", xs[0], spline.SecondDerivativeOn(1, xs[0]))
	fmt.Printf("S''(x_%d = %.2f) = 2c_%d + 6d_%d*h_%d = %.10f\n", n, xs[n], n, n, n, spline.SecondDerivativeOn(n, xs[n]))
	fmt.Println()

	figure := plot.Figure{
		Title:  "Кубический сплайн (естественные граничные условия)",
		XLabel: "x",
		YLabel: "y",
	}
	for i := 1; i <= n; i++ {
		segmentIndex := i
		sx, sy := plot.Sample(func(x float64) float64 { return spline.ValueOnSegment(segmentIndex, x) }, xs[i-1], xs[i], 100)
		figure.Series = append(figure.Series, plot.Series{
			Name: fmt.Sprintf("S_%d на [%.0f, %.0f]", i, xs[i-1], xs[i]), X: sx, Y: sy, Color: plot.Palette[(i-1)%len(plot.Palette)], Width: 2.5,
		})
	}
	figure.Series = append(figure.Series,
		plot.Series{Name: "Узлы (x_i, f_i)", X: xs, Y: fs, Points: true, Color: "#222222"},
		plot.Series{Name: fmt.Sprintf("(X*, S(X*)) = (%.2f, %.4f)", xStar, value), X: []float64{xStar}, Y: []float64{value}, Points: true, Color: "#e377c2", Width: 7},
	)
	plotPath := filepath.Join(filepath.Dir(inputPath), "plot.svg")
	if err := figure.Save(plotPath); err != nil {
		return fmt.Errorf("save plot: %w", err)
	}
	return nil
}
