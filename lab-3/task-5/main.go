package task_5

import (
	"fmt"
	"math"
	"numerical-methods/utils"
)

func Function(x float64) float64 {
	return x / (2*x + 5)
}

func Antiderivative(x float64) float64 {
	return x/2 - 1.25*math.Log(math.Abs(2*x+5))
}

type Method struct {
	Name  string
	Order int
	Apply func(f func(float64) float64, a, b, h float64) (float64, error)
}

func steps(a, b, h float64) (int, error) {
	n := (b - a) / h
	rounded := math.Round(n)
	if h <= 0 || rounded < 1 || math.Abs(n-rounded) > utils.Eps {
		return 0, fmt.Errorf("step h = %g does not divide [%g, %g]", h, a, b)
	}
	return int(rounded), nil
}

func LeftRectangles(f func(float64) float64, a, b, h float64) (float64, error) {
	n, err := steps(a, b, h)
	if err != nil {
		return 0, err
	}
	sum := 0.0
	for i := range n {
		sum += f(a + float64(i)*h)
	}
	return h * sum, nil
}

func RightRectangles(f func(float64) float64, a, b, h float64) (float64, error) {
	n, err := steps(a, b, h)
	if err != nil {
		return 0, err
	}
	sum := 0.0
	for i := 1; i <= n; i++ {
		sum += f(a + float64(i)*h)
	}
	return h * sum, nil
}

func MiddleRectangles(f func(float64) float64, a, b, h float64) (float64, error) {
	n, err := steps(a, b, h)
	if err != nil {
		return 0, err
	}
	sum := 0.0
	for i := range n {
		sum += f(a + (float64(i)+0.5)*h)
	}
	return h * sum, nil
}

func Trapezoids(f func(float64) float64, a, b, h float64) (float64, error) {
	n, err := steps(a, b, h)
	if err != nil {
		return 0, err
	}
	sum := (f(a) + f(b)) / 2
	for i := 1; i < n; i++ {
		sum += f(a + float64(i)*h)
	}
	return h * sum, nil
}

func Simpson(f func(float64) float64, a, b, h float64) (float64, error) {
	n, err := steps(a, b, h)
	if err != nil {
		return 0, err
	}
	if n%2 != 0 {
		return 0, fmt.Errorf("the Simpson method needs an even number of segments, got %d for h = %g", n, h)
	}
	sum := f(a) + f(b)
	for i := 1; i < n; i++ {
		if i%2 == 1 {
			sum += 4 * f(a+float64(i)*h)
		} else {
			sum += 2 * f(a+float64(i)*h)
		}
	}
	return h / 3 * sum, nil
}

func RungeRomberg(f1, f2, h1, h2 float64, p int) (float64, float64) {
	k := h1 / h2
	correction := (f2 - f1) / (math.Pow(k, float64(p)) - 1)
	return f2 + correction, math.Abs(correction)
}

func Run(inputPath string) error {
	data, err := utils.ReadRows(inputPath, 4)
	if err != nil {
		return err
	}
	a, b, h1, h2 := data[0][0], data[0][1], data[0][2], data[0][3]
	if b <= a || h1 <= 0 || h2 <= 0 || h1 == h2 {
		return fmt.Errorf("incorrect integration parameters")
	}
	if 2*a+5 <= 0 && 2*b+5 >= 0 {
		return fmt.Errorf("function has a singularity at x = -2.5 inside the interval")
	}

	methods := []Method{
		{"Левые прямоугольники", 1, LeftRectangles},
		{"Правые прямоугольники", 1, RightRectangles},
		{"Средние прямоугольники", 2, MiddleRectangles},
		{"Трапеции", 2, Trapezoids},
		{"Симпсон", 4, Simpson},
	}
	exact := Antiderivative(b) - Antiderivative(a)

	fmt.Printf("Пределы: X0 = %g, Xk = %g; шаги: h1 = %g, h2 = %g\n", a, b, h1, h2)
	fmt.Printf("Точное значение: F = [x/2 - (5/4) ln|2x+5|] = 1 - 1.25*ln(7/3) = %.10f\n", exact)
	fmt.Println()

	header := []string{"F(h1)", "F(h2)", "p", "F_РР", "delta(h2) по РР", "|F(h1)-F|", "|F(h2)-F|", "|F_РР-F|"}
	widths := []int{14, 14, 2, 14, 12, 11, 11, 11}
	fmt.Print(utils.PadRight("Метод", 23))
	for i, title := range header {
		fmt.Print(" " + utils.PadLeft(title, widths[i]))
	}
	fmt.Println()
	bestName, bestError := "", math.Inf(1)
	bestRRName, bestRRError := "", math.Inf(1)
	for _, method := range methods {
		f1, err := method.Apply(Function, a, b, h1)
		if err != nil {
			return fmt.Errorf("%s: %w", method.Name, err)
		}
		f2, err := method.Apply(Function, a, b, h2)
		if err != nil {
			return fmt.Errorf("%s: %w", method.Name, err)
		}
		refined, rrError := RungeRomberg(f1, f2, h1, h2, method.Order)
		e1, e2, eRR := math.Abs(f1-exact), math.Abs(f2-exact), math.Abs(refined-exact)
		fmt.Printf("%s %14.10f %14.10f %2d %14.10f %12.4e %11.4e %11.4e %11.4e\n", utils.PadRight(method.Name, 23), f1, f2, method.Order, refined, rrError, e1, e2, eRR)
		if e2 < bestError {
			bestName, bestError = method.Name, e2
		}
		if eRR < bestRRError {
			bestRRName, bestRRError = method.Name, eRR
		}
	}
	fmt.Println()
	fmt.Printf("Наиболее точный результат с шагом h2: %s (погрешность %.4e)\n", bestName, bestError)
	fmt.Printf("Наиболее точное уточнённое значение (Рунге–Ромберг): %s (погрешность %.4e)\n", bestRRName, bestRRError)
	return nil
}
