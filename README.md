# Genetic Programming for Electricity Load Prediction

## Introduction

This project implements **Structural-Based Grammatical Evolution (SBGE)**.

The system evolves mathematical expressions that predict future electricity load values from previously observed readings. The model performs symbolic regression, where the algorithm searches for mathematical expressions that minimise prediction error.

The fitness of each evolved individual is evaluated using **Mean Squared Error (MSE)**, while **Mean Absolute Error (MAE)** is also reported for interpretation of prediction accuracy.

---

# Requirements

Before running the program, ensure the following are installed:

- **Java Development Kit (JDK) 11 or newer**
- At least **512 MB available memory**
- A terminal environment (Command Prompt, PowerShell, Bash, or Linux terminal)

Verify installation:


java -version
javac -version


---

# Project Structure

text
Assignment 1/

    ExpressionEvaluator.java
    EnergyPredictSBGE.java
    GEIndividual.java
    BNFGrammar.java
    Residential_Energy_Dataset_UK-2014-2020.csv
    Makefile
    run.bat
    README.md


> The dataset file must remain in the same directory as the Java source files.

---

# How to Run (Linux / macOS)

## Compile All Files


make compile


## Run Structural-Based Grammatical Evolution (SBGE)


make run-sbge


## Clean Compiled Files


make clean


---

# How to Run (Windows)

## Using run.bat (Recommended)


.\run.bat

Select mode:
1. Structural-Based GE (SBGE)
Enter choice (1): 

---
1
## Manual Compilation


javac ExpressionEvaluator.java EnergyPredictSBGE.java GEIndividual.java BNFGrammar.java


## Run SBGE


java EnergyPredictSBGE


---

# Parameter Configuration

All algorithm parameters can be modified directly inside:

text
EnergyPredictSBGE.java


## Main Parameters

Parameter | Description |

`POPSIZE` = Population size 
`MAXGEN` = Number of generations 
`CODON_LENGTH` = Number of codons per individual 
`MUTATION_LOW` = Low mutation probability 
`MUTATION_HIGH` = High mutation probability 
`TOURNAMENT_SIZE` = Tournament selection size 
`SUCCESS_THRESHOLD` = Target fitness threshold 

---

## Dataset Size Control (IMPORTANT)

The parameter controlling how much data is loaded is:

java
static int numForDataset = 0;


### Dataset Modes

| Value | Description |
 `0` = Loads half of the dataset 
 `1` = Loads 10,000 rows 
 `2+` = Loads the full dataset 


### Performance and Runtime Considerations

The size of the dataset directly impacts execution time:

- **Larger datasets** = more accurate results but slower execution
- **Smaller datasets** = faster execution for testing and debugging

**Recommendation:**
- Use `numForDataset = 0` or `1` during development and testing
- Use `numForDataset >= 2` for final evaluation
### Purpose of Dataset Control

The `numForDataset` parameter was introduced to address long runtime issues. By allowing partial dataset loading, the system enables faster experimentation while maintaining the ability to scale to the full dataset when required.



## If No Improvement Is Observed

If the algorithm does not improve over generations, consider adjusting the following parameters:

- Mutation rate (`MUTATION`)
- Crossover rate (`CROSSOVER`)
- Population size (`PopSize`)
- Maximum tree depth (`MaxDepth`)
- Tournament size (`TOURNAMENT`)
- Dataset size (`numForDataset`)

Reducing the dataset size can help to:

- Identify issues faster
- Improve debugging speed
- Reduce long execution times

---

# Dataset Information

The dataset used is:

text
Residential_Energy_Dataset_UK-2014-2020.csv


Only the `Electricity_load` column is used for prediction.

The system creates time-series windows of previous electricity load values to predict the next load value.

---

# Adaptive Mutation and Crossover

The system uses adaptive mutation and crossover probabilities based on population similarity.

- When population similarity is low, mutation remains low to preserve diversity.
- When population similarity becomes high, mutation increases to introduce new variations and avoid convergence.

Crossover probability is calculated from the adaptive mutation probability.

This helps balance:
- Exploration of new solutions
- Exploitation of good existing solutions

---

# Expected Output

During execution, the program prints:

- Generation progress
- Best fitness values
- Population statistics
- Prediction results
- Final evaluation summary

---

# Notes

- Larger datasets may increase runtime significantly.
- Mutation and crossover probabilities can strongly affect performance.
- The dataset file must remain in the same folder as the Java source files.
- The system preserves time order in the dataset to avoid look-ahead leakage during prediction.

---
