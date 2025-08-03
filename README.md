# simcol
## Program for Simulating Collusion
Simulate Collusion based on Simulation of industries with collusive and non-collusive periods over time, detected and undetected, based on Harrington/Chang (2015). This simulation is described in detail in the paper "Simulating Collusion and Enforcement" by Bellert and Günster (2025).

## Requirements
Install R: https://cran.rstudio.com/

Install RStudio: https://posit.co/downloads/

In your RStudio console:

`library(devtools)`

`install_github("belln1/simcol");`

`install.packages(c("dplyr", "EnvStats"))`


## Getting Started
Type `?simcol` and run the described examples.

Run R Markdown with code examples for simulation, summary statistics, plot and estimates:  
examples/examples.Rmd

Run R file simulating different Models 1, 2, 2B, 3:  
examples/test_simulate_collusion.R


## Questions and comments 
bell@zhaw.ch
