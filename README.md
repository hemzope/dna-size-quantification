# DNA Size Quantification and Suspect Matching

A Python-based DNA fragment size quantification tool using gel electrophoresis data.

This software creates a DNA ladder calibration curve, estimates fragment sizes of unknown samples, generates plots and CSV reports, and automatically ranks suspects based on similarity to an evidence sample.

---

## Features

* DNA ladder calibration using linear regression
* Automatic calculation of:

  * Slope
  * Intercept
  * R² (coefficient of determination)
* DNA fragment size estimation
* Automatic suspect matching
* CSV export of all results
* PNG export of all plots
* Supports any number of suspects
* JSON-based input file
* Automatic ranking of suspects by similarity to evidence

---

## Folder Structure

```text
project/
│
├── main.py
├── data.json
├── sample_data.json
├── requirements.txt
├── LICENSE
├── README.md
├── .gitignore
│
├── sample_results/
│   └── ...
│
└── results/
```

The `results` folder is created automatically when the program runs.

After execution:

```text
results/
│
├── ladder/
│   ├── ladder_data.csv
│   └── standard_curve.png
│
├── evidence/
│   ├── evidence_results.csv
│   └── evidence_plot.png
│
├── suspects/
│   ├── suspect1_results.csv
│   ├── suspect1_plot.png
│   ├── suspect2_results.csv
│   ├── suspect2_plot.png
│   └── ...
│
└── comparison/
    ├── suspect_ranking.csv
    └── match_report.csv
```

---

## Requirements

### Python Version

This project is tested on `Python 3.14`.
However, it is expected to work on `Python 3.10+`.

Check your Python version:

```bash
python --version
```

or

```bash
python3 --version
```

---

## Installation

### 1. Clone the Repository

```bash
git clone <repository-url>
cd <repository-name>
```

### 2. Install Dependencies

A requirements file is included.

Install all required packages:

```bash
pip install -r requirements.txt
```

---

## Input Data

All user-editable data is stored in:

```text
data.json
```

No modification of the Python code is required.

---

## Measurement Guidelines

### Distance Units

Measure all DNA band migration distances using the same unit.

Recommended:

```text
Millimetres (mm)
```

Measure from:

```text
Bottom of loading well
        ↓
Centre of DNA band
```

Use the same measurement method for:

* DNA ladder
* Evidence sample
* All suspects

---

## JSON Format

Example:

```json
{
    "ladder": {
        "bp": [3000, 1500, 1000, 900, 800, 700, 600, 500, 400, 300, 200],
        "distance": [262, 397, 518, 552, 590, 624, 682, 740, 804, 880, 965]
    },

    "evidence": {
        "distance": [365, 465, 585, 828]
    },

    "suspects": {
        "alice": {
            "distance": [272, 373, 430, 730, 833]
        },
        "bob": {
            "distance": [335, 445, 573, 822]
        },
        "charlie": {
            "distance": [340, 445, 570, 610, 827]
        }
    }
}
```

---

## Running the Program

Execute:

```bash
python main.py
```

or

```bash
python3 main.py
```

The program will automatically:

1. Read the JSON file
2. Create the ladder calibration curve
3. Estimate DNA fragment sizes
4. Generate plots
5. Generate CSV files
6. Rank suspects
7. Save all results

---

## Understanding the Outputs

### Standard Curve

The DNA ladder is used to generate a calibration curve:

```text
log10(bp) = (m × distance) + b
```

This relationship is then used to estimate fragment sizes of unknown samples.

---

### R² Value

The coefficient of determination (R²) indicates how well the ladder data fits the calibration line.

Values closer to:

```text
1.000
```

indicate a better fit.

---

### Suspect Ranking

Suspects are automatically ranked according to similarity to the evidence sample.

Example:

```text
Suspect              Error
--------------------------------
Alice                0.0021
Bob                  0.0548
Charlie              0.1137
```

Lower values indicate a closer match.

---

## Included Files

### Sample Data

A sample input file is included:

```text
sample_data.json
```

This can be used to test the software before entering your own measurements.

### Sample Results

Example outputs are included in:

```text
sample_results/
```

These demonstrate the expected CSV files and plots.

---

## Scientific Notes

This software assumes:

* A linear relationship between migration distance and log10(fragment size).
* Measurements were collected consistently.
* Ladder and unknown samples were measured using the same method.

This project is intended primarily for:

* Educational demonstrations
* Science workshops
* Student projects
* Introductory forensic analysis exercises

It should not be used for legal, forensic, clinical, or research decision-making without independent validation.

---

## License

This project is distributed under the MIT License, included in the repository.

Please read the `LICENSE` file before using, modifying, or redistributing this software.

By using this software, you agree to comply with the terms of that license.

---

## Acknowledgements

This project was developed as an educational tool to demonstrate the principles of:

* Gel electrophoresis
* DNA fragment size estimation
* Regression analysis
* Scientific data processing using Python
