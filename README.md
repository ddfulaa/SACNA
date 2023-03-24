<p align="center">
  <img width="290" height="225" src="https://github.com/ddfulaa/SACNA/blob/main/sacna_logo.png?raw=true">
</p>

# SACNA (Semi-Algebraic Chemical Network Analyzer)

SACNA is a Mathematica package developed by Daniel Fula and J. Montoya for the analysis and simulation of chemical reaction networks that exhibit chiral amplifiers, using the Collins' algorithm. The package provides 2 ways of analysis and a simulator.

## Requirements

* Mathematica 12

## How to use it
The Tutorial folder includes several Mathematica® notebooks that illustrate the usage of SACNA in different chiral chemical reaction networks. These examples represent possible results when using SACNA to analyze chiral networks:

Calvin Model: The algorithm terminates without finding a solution, determining that the network has no chiral amplifier. (Calvin model)
Calvin-LES Model: The algorithm finds solutions for all Routh-Hurwitz conditions. (Calvin-LES model)
Hochberg-Ribo Model: The algorithm finds solutions for one of the Routh-Hurwitz conditions. (Hochberg-Ribo replicator model)
APED Model: The entered network is too complex, making most of the Routh-Hurwitz conditions impossible to analyze or have no short-term solutions. (APED model)
APED-reversible Model: The entered network is too complex, making all Routh-Hurwitz conditions impossible to analyze or have no short-term solutions. (APED reversible model)
Please note that the examples are located in the Tutorial folder and are accompanied by their respective PDF files with explanations of the package usage.

## License
This code is released under Apache 2.0
