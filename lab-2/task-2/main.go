package task_2

import (
	"fmt"
	"math"
	"numerical-methods/utils"
)

const maxIterations = math.MaxInt32

type Point struct {
	x1 float64
	x2 float64
}

type Iteration struct {
	Number        int
	Point         Point
	deltaX        float64
	ErrorEstimate float64
}

// region Helpers
func printHistory(history []Iteration) {
	fmt.Printf("%3s  %17s  %17s  %17s  %17s\n", "k", "x1^(k)", "x2^(k)", "max |delta x|", "Error estimate")
	for _, item := range history {
		fmt.Printf("%3d  %17.10f  %17.10f  %17.10e  %17.10e\n", item.Number, item.Point.x1, item.Point.x2, item.deltaX, item.ErrorEstimate)
	}
}

//endregion Helpers

func Function1(point Point) float64 {
	return (point.x1*point.x1+4)*point.x2 - 8
}

func Function2(point Point) float64 {
	return (point.x1-1)*(point.x1-1) + (point.x2-1)*(point.x2-1) - 4
}

func MaxDifference(first, second Point) float64 {
	return math.Max(math.Abs(first.x1-second.x1), math.Abs(first.x2-second.x2))
}

func IterationFunction(point Point) (Point, error) {
	f := 4 - (point.x2-1)*(point.x2-1)
	if f < 0 {
		return Point{}, fmt.Errorf("iteration function is undefined for x2 = %.8f", point.x2)
	}

	denominator := point.x1*point.x1 + 4
	if denominator == 0 {
		return Point{}, fmt.Errorf("iteration function denominator is zero")
	}

	return Point{x1: 1 + math.Sqrt(f), x2: 8 / denominator}, nil
}

func PhiDerivative() float64 {
	phi1Derivative := 0.5 / math.Sqrt(4-0.25)
	phi2Derivative := 16 * 2.8 / math.Pow(2.8*2.8+4, 2)
	return math.Max(phi1Derivative, phi2Derivative)
}

func SimpleIterations(eps float64) (Point, []Iteration, error) {
	if eps <= 0 {
		return Point{}, nil, fmt.Errorf("eps must be positive")
	}

	prev := Point{x1: 3, x2: 0.6}
	q := PhiDerivative()
	history := make([]Iteration, 0)

	for k := 1; k <= maxIterations; k++ {
		cur, err := IterationFunction(prev)
		if err != nil {
			return Point{}, nil, err
		}

		deltaX := MaxDifference(cur, prev)
		errorEstimate := q / (1 - q) * deltaX
		history = append(history, Iteration{
			Number:        k,
			Point:         cur,
			deltaX:        deltaX,
			ErrorEstimate: errorEstimate,
		})

		if errorEstimate <= eps {
			return cur, history, nil
		}

		prev = cur
	}

	return Point{}, nil, fmt.Errorf("simple iterations didt finish after %d iterations", maxIterations)
}

func Newton(eps float64) (Point, []Iteration, error) {
	if eps <= 0 {
		return Point{}, nil, fmt.Errorf("eps must be positive")
	}

	cur := Point{x1: 3, x2: 0.6}
	history := make([]Iteration, 0)

	for k := 1; k <= maxIterations; k++ {
		f1 := Function1(cur)
		f2 := Function2(cur)
		j11 := 2 * cur.x1 * cur.x2
		j12 := cur.x1*cur.x1 + 4
		j21 := 2 * (cur.x1 - 1)
		j22 := 2 * (cur.x2 - 1)
		determinant := j11*j22 - j12*j21
		if determinant == 0 {
			return Point{}, nil, fmt.Errorf("the Jacobian determinant is zero at (%.8f, %.8f)", cur.x1, cur.x2)
		}

		deltaX1 := (-f1*j22 + f2*j12) / determinant
		deltaX2 := (j21*f1 - j11*f2) / determinant
		next := Point{x1: cur.x1 + deltaX1, x2: cur.x2 + deltaX2}
		deltaX := MaxDifference(next, cur)
		history = append(history, Iteration{
			Number:        k,
			Point:         next,
			deltaX:        deltaX,
			ErrorEstimate: math.NaN(),
		})

		if deltaX <= eps {
			return next, history, nil
		}

		cur = next
	}

	return Point{}, nil, fmt.Errorf("the Newton method didt finish after %d iterations", maxIterations)
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

	fmt.Println("Task:")
	fmt.Println("(x1^2 + 4) * x2 - 8 = 0")
	fmt.Println("(x1 - 1)^2 + (x2 - 1)^2 - 4 = 0")
	fmt.Println("Positive solution: x1 in [2.8, 3.1], x2 in [0.5, 0.75]")
	fmt.Println("Epsilon:", epsilon)
	fmt.Println()
	fmt.Println("SimpleIterationsMethod:")
	fmt.Println("x^(0) = (3, 0.6)")
	fmt.Println("phi1(x2) = 1 + sqrt(4 - (x2 - 1)^2)")
	fmt.Println("phi2(x1) = 8 / (x1^2 + 4)")
	fmt.Printf("q = %.10f\n", PhiDerivative())
	printHistory(simpleIterationsHistory)
	fmt.Printf("Solution: (%.10f, %.10f)\n", simpleIterationsSolution.x1, simpleIterationsSolution.x2)
	fmt.Printf("Iterations: %d\n", len(simpleIterationsHistory))
	fmt.Println()
	fmt.Println("NewtonMethod:")
	fmt.Println("x^(0) = (3, 0.6)")
	printHistory(newtonMethodHistory)
	fmt.Printf("Solution: (%.10f, %.10f)\n", newtonMethodSolution.x1, newtonMethodSolution.x2)
	fmt.Printf("Iterations: %d\n", len(newtonMethodHistory))

	return nil
}
