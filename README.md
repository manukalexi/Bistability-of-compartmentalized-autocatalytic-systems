Repository for the software described in

# Bistability of compartmentalized autocatalytic systems at the dawn of hereditary variation

Alexa Iván<sup>1,2,4</sup> ; András Szilágyi<sup>1,2</sup> ; Eörs Szathmáry<sup>1,2,3</sup>
 
<sup>1</sup> Center for the Conceptual Foundations of Science, Parmenides Foundation, Pöcking, Germany   
<sup>2</sup> Institute of Evolution, HUN-REN Centre for Ecological Research, Budapest, Hungary   
<sup>3</sup> Department of Plant Systematics, Ecology and Theoretical Biology, ELTE Eötvös Loránd University, Budapest, Hungary   
<sup>4</sup> Doctoral School of Biology, ELTE Eötvös Loránd University, Budapest, Hungary   

## Models

This repository contains the models used to investigate alternative stable chemical states in interacting autocatalytic networks, both without and within compartments. The models describe two mutually interacting autocatalytic cycles, A and B, producing molecules X and Y, respectively. Their coupled dynamics can give rise to bistability and hysteresis, with the resulting state depending on the initial conditions.

The chemostat (CSTR) model describes the reaction dynamics in a continuously supplied and diluted environment, while the compartment-level selection model implements the same dynamics within a population of reproducing compartments. Molecules are stochastically partitioned at division, generating variation in compartment composition. Environmental change is represented by two environments, `E_X` and `E_Y`, which favor compartments enriched in X or Y, respectively, through composition-dependent autocatalytic rates. This can drive regime shifts either through within-compartment dynamics or, when the system remains within the hysteresis region, through compartment-level selection acting on stochastically generated variants.

The repository also includes extended reaction networks with explicitly represented inhibitory molecules, modified fragmentation rules, and longer autocatalytic cycles, under both reversible and irreversible reaction schemes and periodically changing environments.

## Repository structure

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
./Fig.2
```
```bash
./Fig.3
```
```bash
./Fig.S1
```
```bash
./Fig.S2
```
```bash
./Fig.S3S4
```
```bash
./Fig.S5
```

4. Save the program output to a file using output redirection (`>`):

```bash
./Fig.2 > output_filename.txt
```
```bash
./Fig.3 > output_filename.txt
```
```bash
./Fig.S1 > output_filename.txt
```
```bash
./Fig.S2 > output_filename.txt
```
```bash
./Fig.S3S4 > output_filename.txt
```
```bash
./Fig.S5 > output_filename.txt
```
5. Use the raw output files for further data evaluation and visualization.

6. For the Mathematica codes, open the corresponding files in Mathematica and evaluate the code directly; no compilation is required.

## Parameters of the models

### General parameters

- `a`: autocatalytic rate constant of cycle A
- `b`: autocatalytic rate constant of cycle B
- `g`: annihilation rate constant for the destructive reaction between molecules X and Y
- `c`: parameter controlling the strength of the composition-dependent environmental effect
- `m`: interconversion rate constant between molecules X and Y
- `x`: number of X molecules
- `y`: number of Y molecules
- `S`: total number of molecules required for compartment division
- `N`: total number of compartments in the population
- `D`: number of daughter compartments produced upon division
- `T`: number of additional population members replaced by daughter compartments during division (one daughter always replaces the parent compartment, while the remaining `T` daughters replace randomly selected compartments)
- `τ`: period of environmental change

### Fig. 1

- `r`: fixed concentration of the resource R
- `f`: outflow-like competition parameter preventing unbounded autocatalyst growth

### Fig. 2

**Chemostat (CSTR) model:**

- `d`: dilution rate constant of the CSTR
- `R`: resource concentration in the inflow
- `φ`: outflow/turnover parameter controlling the chemostat dilution

**Compartmentalized model:**

- `R_0`: initial resource concentration

### Fig. 3

**Chemostat (CSTR) model:**

- `R_0`: initial resource concentration
- `d`: dilution rate constant of the chemostat

**Compartmentalized model:**

- `R`: resource concentration
- `m_x`: interconversion rate constant from molecule X to Y
- `m_y`: interconversion rate constant from molecule Y to X

### Fig. S1

- `g_x`, `g_y`: reaction rate constants for the formation of inert complexes between molecules X and W, and Y and V, respectively
- `h_x`, `h_y`: conversion rate constants of molecules X and Y into the inhibitory molecules V and W, respectively
- `a'`, `b'`, `h_x'`, `h_y'`, `g_x'`, `g_y'`, `m'`: reverse reaction rate constants, set to 10% of the corresponding forward rates

### Fig. S2, S3, S4

- `g_x`, `g_y`: reaction rate constants for the formation of inert complexes between molecules X and W, and Y and V, respectively
- `h_x`, `h_y`: conversion rate constants of molecules X and Y into the inhibitory molecules V and W, respectively

### Fig. S5

- `a_1`, `a_2`: autocatalytic rate constants of the two steps of cycle A
- `b_1`, `b_2`: autocatalytic rate constants of the two steps of cycle B
- `g_x`, `g_y`: reaction rate constants for the formation of inert complexes between molecules X and W, and Y and V, respectively
- `h_x`, `h_y`: conversion rate constants of molecules X and Y into the inhibitory molecules V and W, respectively

## Disclaimer
**The code provided in this repository is released for research purposes only and is provided without warranty or liability.**
