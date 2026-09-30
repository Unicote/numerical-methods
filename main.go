package main

import (
	"errors"
	"flag"
	"fmt"
	"numerical-methods/lab-1/task-1"
	"numerical-methods/lab-1/task-2"
	"numerical-methods/lab-1/task-3"
	"numerical-methods/lab-1/task-4"
	"numerical-methods/lab-1/task-5"
	"os"
	"path/filepath"
)

type taskRunner func(string) error

func findTask(lab, task int) (taskRunner, error) {
	if lab != 1 {
		return nil, fmt.Errorf("laboratory work %d is not implemented", lab)
	}

	switch task {
	case 1:
		return task_1.Run, nil
	case 2:
		return task_2.Run, nil
	case 3:
		return task_3.Run, nil
	case 4:
		return task_4.Run, nil
	case 5:
		return task_5.Run, nil
	default:
		return nil, fmt.Errorf("task %d of laboratory work %d is not implemented", task, lab)
	}
}

func main() {
	lab := flag.Int("lab", 0, "laboratory work number")
	task := flag.Int("task", 0, "task number")
	input := flag.String("input", "", "path to input file")
	flag.Parse()

	if *lab <= 0 || *task <= 0 {
		_, err := fmt.Fprintln(os.Stderr, "Usage: go run . -lab 1 -task 1")
		if err != nil {
			return
		}
		flag.PrintDefaults()
		os.Exit(2)
	}

	run, err := findTask(*lab, *task)
	if err != nil {
		_, err := fmt.Fprintln(os.Stderr, "Error:", err)
		if err != nil {
			return
		}
		os.Exit(1)
	}

	inputPath := filepath.Join(fmt.Sprintf("lab-%d", *lab), fmt.Sprintf("task-%d", *task), "input.txt")
	if *input != "" {
		inputPath = *input
	}

	if err := run(inputPath); err != nil {
		if errors.Is(err, os.ErrNotExist) {
			_, err := fmt.Fprintf(os.Stderr, "Error: input file %q was not found\n", inputPath)
			if err != nil {
				return
			}
		} else {
			_, err := fmt.Fprintln(os.Stderr, "Error:", err)
			if err != nil {
				return
			}
		}
		os.Exit(1)
	}
}
