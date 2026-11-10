package plot

import (
	"bytes"
	_ "embed"
	"encoding/json"
	"fmt"
	"os/exec"
)

//go:embed plot.py
var script string

type Series struct {
	Name   string    `json:"name"`
	X      []float64 `json:"x"`
	Y      []float64 `json:"y"`
	Points bool      `json:"points"` // точки вместо линии
	Color  string    `json:"color"`
	Dashed bool      `json:"dashed"`
	Width  float64   `json:"width"`
}

type Figure struct {
	Title  string
	XLabel string
	YLabel string
	Series []Series
}

var Palette = []string{"#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf"}

func Sample(f func(float64) float64, a, b float64, n int) ([]float64, []float64) {
	xs := make([]float64, n+1)
	ys := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + (b-a)*float64(i)/float64(n)
		ys[i] = f(xs[i])
	}
	return xs, ys
}

type request struct {
	Path   string   `json:"path"`
	Title  string   `json:"title"`
	XLabel string   `json:"xlabel"`
	YLabel string   `json:"ylabel"`
	Series []Series `json:"series"`
}

func (f Figure) Save(path string) error {
	data := request{
		Path:   path,
		Title:  f.Title,
		XLabel: f.XLabel,
		YLabel: f.YLabel,
		Series: f.Series,
	}

	payload, err := json.Marshal(data)
	if err != nil {
		return fmt.Errorf("encode plot data: %w", err)
	}

	cmd := exec.Command("python3", "-c", script)
	cmd.Stdin = bytes.NewReader(payload)
	if output, err := cmd.CombinedOutput(); err != nil {
		return fmt.Errorf("python3 error: %w\n%s", err, output)
	}
	return nil
}
