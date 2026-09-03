Repository for the software described in

# Bistability of compartmentalized autocatalytic systems at the dawn of hereditary variation

Alexa Iván<sup>1,2,4</sup>; András Szilágyi<sup>1,2</sup>; Eörs Szathmáry<sup>1,2,3</sup>
 
<sup>1</sup> Center for the Conceptual Foundations of Science, Parmenides Foundation, Pöcking, Germany   
<sup>2</sup> Institute of Evolution, HUN-REN Centre for Ecological Research, Budapest, Hungary   
<sup>3</sup> Department of Plant Systematics, Ecology and Theoretical Biology, ELTE Eötvös Loránd University, Budapest, Hungary   
<sup>4</sup> Doctoral School of Biology, ELTE Eötvös Loránd University, Budapest, Hungary   

## Repository Structure

The repository is organized according to the figures associated with the computational analyses:

- `Fig.1/` — Wolfram Language code
- `Fig.2/` — C code
- `Fig.3/` — Wolfram Language and C code
- `Fig.S1/` — C code
- `Fig.S2/` — C code
- `Fig.S3S4/` — C code
- `Fig.S5/` — C code

## Installation and Usage

1. Download the directory corresponding to the figure of interest. Each directory contains the complete C source code and/or Wolfram Language code required to reproduce the associated computational analysis.

2. On Linux/Ubuntu systems with the GCC toolchain installed, navigate to the directory or subdirectory of the selected figure and compile the C code using:

```bash
make
```

3. Run the compiled executable corresponding to the selected model:

```bash
./Fig.2.c
```
```bash
./Fig.S1.c
```
```bash
./Fig.S2.c
```
```bash
./Fig.S3S4.c
```
```bash
./Fig.S5.c
```

4. Save the program output to a file using output redirection (`>`):

```bash
./Fig.2.c > output_filename.txt
```
```bash
./Fig.S1.c > output_filename.txt
```
```bash
./Fig.S2.c > output_filename.txt
```
```bash
./Fig.S3S4.c > output_filename.txt
```
```bash
./Fig.S5.c > output_filename.txt
```
5. Use the raw output files for further data evaluation and visualization.

6. For the Mathematica codes, open the corresponding files in Mathematica and evaluate the code directly; no compilation is required.

## Parameters of the model

### Model Parameters

### Fig.1

- `a`: autocatalytic rate constant of cycle A
- `b`: autocatalytic rate constant of cycle B
- `r`: fixed concentration of the resource R
- `m`: interconversion rate constant between molecules X and Y
- `g`: annihilation rate constant for the destructive reaction between molecules X and Y
- `f`: outflow-like competition parameter preventing unbounded autocatalyst growth
- `x`: concentration of molecule X
- `y`: concentration of molecule Y
- `t`: time

### Fig.2

Common parameters:

- `a`: autocatalytic rate constant of cycle A
- `b`: autocatalytic rate constant of cycle B
- `m`: interconversion rate constant between molecules X and Y
- `g`: annihilation rate constant for the destructive reaction between molecules X and Y

Chemostat (CSTR) model:

- `d`: dilution rate constant of the CSTR
- `R`: resource concentration in the inflow
- `φ`: outflow/turnover parameter controlling the chemostat dilution 

Compartmentalized model:

- `R_0`: initial resource concentration
- `c`: parameter controlling the strength of the composition-dependent environmental effect
- `S`: total number of molecules required for compartment division, e.g. `x[i] + y[i] + p[i] = S`
- `N`: total number of compartments in the population
- `D`: number of daughter compartments produced upon division
- `T`: number of additional population members replaced by daughter compartments during division (one daughter always replaces the parent compartment, while the remaining `T` daughters replace randomly selected compartments)

### Fig.3

Chemostat (CSTR) model:

- `R_0`: initial resource concentration
- `a`: autocatalytic rate constant of cycle A
- `b`: autocatalytic rate constant of cycle B
- `m`: interconversion rate constant between molecules X and Y
- `g`: annihilation rate constant for the destructive reaction between molecules X and Y
- `d`: dilution rate constant of the chemostat
- `c`: parameter controlling the strength of the composition-dependent environmental effect

Compartmentalized model:

- `R`: resource concentration
- `a`: autocatalytic rate constant of cycle A
- `b`: autocatalytic rate constant of cycle B
- `m_x`: interconversion rate constant from molecule X to Y
- `m_y`: interconversion rate constant from molecule Y to X
- `g`: annihilation rate constant for the destructive reaction between molecules X and Y
- `c`: parameter controlling the strength of the composition-dependent environmental effect
- `S`: critical total number of molecules required for compartment division
- `N`: total number of compartments in the population
- `D`: number of daughter compartments produced upon division
- `T`: number of additional population members replaced by daughter compartments during division (one daughter always replaces the parent compartment, while the remaining `T` daughters replace randomly selected compartments)
 
## Disclaimer
**The code provided in this repository is released for research purposes only and is provided without warranty or liability.**
