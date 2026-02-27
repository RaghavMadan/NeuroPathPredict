# Contributing to NeuroPathPredict Pipeline V1.1

Thank you for your interest in contributing. Please follow these guidelines.

## How to contribute

- **Bug reports and feature requests:** Open a GitHub Issue.
- **Code changes:** Fork the repository, create a branch, make your changes, then open a Pull Request.

## Development setup

1. Clone the repository and open the project root (the directory containing `config/`, `scripts/`, `input/`, `output/`).
2. Copy `config/input_variables.example.txt` to `config/input_variables.txt` and set paths to your data (see README).
3. Copy `config/subject_list.example.txt` to `config/subject_list.txt` and add your subject IDs, or use the example list.
4. Install R dependencies (see README). Run `Rscript scripts/setup.R` to install required packages.
5. Run the pipeline: `Rscript scripts/run_pipeline.R` from the project root.

## Code style

- Use consistent formatting in R and Python scripts.
- Comment non-obvious steps. Keep paths relative to the project root where possible.

## Pull requests

- Describe what changed and why.
- Ensure the pipeline still runs from the project root with the example config (or document any new required setup).

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
