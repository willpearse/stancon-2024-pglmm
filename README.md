# StanCon 2024: Biodiversity Workshop Materials

Welcome to the materials repository for the **StanCon 2024: Biodiversity Workshop**! This repository contains all the resources you'll need to participate in the workshop, including the presentation slides, a manuscript detailing our work, and an R script for hands-on exercises.

## Table of Contents

- [Presentation](#presentation)
- [Manuscript](#manuscript)
- [R Script](#r-script)
- [Getting Started](#getting-started)
- [License](#license)

## Presentation

You can access the presentation slides for the workshop [here](https://docs.google.com/presentation/d/1rT1qQhGu_VumPHMklpg06xn8XCAwQV-sVGNdm69U3h0/edit?usp=sharing).

## Manuscript

This workshop covers an approach that is outlined in detail in [Gallinat & Pearse (2021; Oikos 130(5): 669-79)](https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/oik.08048); you can either follow that link or the PDF is included in this repository ([`pglmm-paper.pdf`](https://github.com/willpearse/stancon-2024-pglmm/blob/main/pglmm-paper.pdf)).

## R script

The R script we will be working through in the workshop (which contains all the `stan` code you need) is available within this repository (([`pglmm-script.R`](https://github.com/willpearse/stancon-2024-pglmm/blob/main/pglmm-script.R))). It's the supplementary materials for the paper described above.

## Getting Started

To get started with the workshop:

1. Clone this repository to your local machine; you can probably do that from a terminal by typing `git clone https://github.com/willpearse/stancon-2024-pglmm.git` or else download everything manually from here.
2. Install the necessary R dependencies for the course; something like `install.packages("rstan")` (if you haven't already installed `stan`) and `install.packages("pez")` is all you need.
3. Open the `R` script outlined in the section above that's now on your computer. Work your way through that!
4. Follow the steps outlined in the talk and also projected on the wall right now. Those are also repeated below.

## What to do on the day

As is hopefully projected on the board right now:

1. Sit, breathe, digest... then ask questions :D
2. Open up the .R file. ‘Question 1’ maps onto ‘Question 1’ in the attached paper, which also maps onto the talk I just gave
3. Run through the example in the .R file with simulated data.
4. Run the model through with real data:
   a. library(pez)
   b. data(laja)
   c. ...build it from there
5. Question 2 is all about competition among species. Can you fit those models
6. Question 3 is all about co-evolution among groups (host vs. viruses). Can you fit those models?
