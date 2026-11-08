package task_1

import (
	"fmt"
	"math"
	"numerical-methods/utils"
)

const maxIterations = math.MaxInt32

// region Helpers
func printHistory(history []Iteration) {
	fmt.Printf("%3s  %17s  %17s  %17s\n", "k", "x^(k)", "|x^(k)-x^(k-1)|", "Error estimate")
	for _, item := range history {
		fmt.Printf("%3d  %17.10f  %17.10e  %17.10e\n", item.Number, item.X, item.deltaX, item.ErrorEstimate)
	}
}

//endregion Helpers

type Iteration struct {
	Number        int
	X             float64
	deltaX        float64
	ErrorEstimate float64
}

func Function(x float64) float64 {
	return math.Pow(2, x) - x*x - 0.5
}

func FunctionDerivative(x float64) float64 {
	return math.Ln2*math.Pow(2, x) - 2*x
}

func IterationFunction(x float64) (float64, error) {
	f := math.Pow(2, x) - 0.5
	result := math.Sqrt(f)
	if math.IsNaN(result) {
		return 0, fmt.Errorf("iteration function is undefined for x = %.8f", x)
	}
	return result, nil
}

func PhiDerivative() float64 {
	return 2 * math.Ln2 / math.Sqrt(3.5)
}

func SimpleIterations(eps float64) (float64, []Iteration, error) {
	if eps <= 0 {
		return 0, nil, fmt.Errorf("epsilon must be positive")
	}

	const x0 = 1.5
	q := PhiDerivative()
	if q >= 1 {
		return 0, nil, fmt.Errorf("q = %.10f >= 1", q)
	}
	prev := x0
	history := make([]Iteration, 0)

	for k := 1; k <= maxIterations; k++ {
		cur, err := IterationFunction(prev)
		if err != nil {
			return 0, nil, err
		}

		deltaX := math.Abs(cur - prev)
		errorEstimate := q / (1 - q) * deltaX
		history = append(history, Iteration{
			Number:        k,
			X:             cur,
			deltaX:        deltaX,
			ErrorEstimate: errorEstimate,
		})

		if errorEstimate <= eps {
			return cur, history, nil
		}

		prev = cur
	}

	return 0, nil, fmt.Errorf("simple iterations didt finish after %d iterations", maxIterations)
}

func Newton(eps float64) (float64, []Iteration, error) {
	if eps <= 0 {
		return 0, nil, fmt.Errorf("eps must be positive")
	}

	const x0 = 1.7
	cur := x0
	history := make([]Iteration, 0)

	for k := 1; k <= maxIterations; k++ {
		derivative := FunctionDerivative(cur)
		if derivative == 0 {
			return 0, nil, fmt.Errorf("derivative is zero when x = %.8f", cur)
		}

		next := cur - Function(cur)/derivative
		deltaX := math.Abs(next - cur)
		history = append(history, Iteration{
			Number:        k,
			X:             next,
			deltaX:        deltaX,
			ErrorEstimate: math.NaN(),
		})

		if deltaX <= eps {
			return next, history, nil
		}

		cur = next
	}

	return 0, nil, fmt.Errorf("the Newton method didt finish after %d iterations", maxIterations)
}

func Run(inputPath string) error {
	epsilon, err := utils.ParseEpsilon(inputPath)
	if err != nil {
		return err
	}

	simpleIterationsSolution, simpleIterationsHistory, err := SimpleIterations(epsilon)
	if err != nil {
		return err
	}

	newtonMethodSolution, newtonMethodHistory, err := Newton(epsilon)
	if err != nil {
		return err
	}

	fmt.Println("Equation: 2^x - x^2 - 0.5 = 0")
	fmt.Println("Positive root: x in [1, 2]")
	fmt.Println("Epsilon:", epsilon)
	fmt.Println()
	fmt.Println("Simple iterations method:")
	fmt.Println("x^(0) = 1.5")
	fmt.Println("phi(x) = sqrt(2^x - 0.5)")
	fmt.Printf("q = %.10f\n", PhiDerivative())
	printHistory(simpleIterationsHistory)
	fmt.Printf("Root: %.10f\n", simpleIterationsSolution)
	fmt.Printf("Iterations : %d\n", len(simpleIterationsHistory))
	fmt.Println()
	fmt.Println("Newton metod:")
	fmt.Println("x^(0) = 1.7")
	printHistory(newtonMethodHistory)
	fmt.Printf("Root: %.10f\n", newtonMethodSolution)
	fmt.Printf("Iterations: %d\n", len(newtonMethodHistory))

	return nil
}
