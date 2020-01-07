# A Framework for Representing Language Acquisition in a Population Setting

Dynamical systems framework (cf Niyogi &amp; Berwick 1997, 2009) for modeling language change in populations of learners extended to represent arbitrary social network structures.

Presented at ACL 2018
https://www.aclweb.org/anthology/P18-1106/

## Setup

Just run **setup.sh** or create a plots directory in the parent directory. Nothing to install. 

## Dependencies

Everything was implemented in Python 2.7. Numpy is the only external library required to run things.

## Experiment Files

There is no user interface yet. Files of the form exp_*.py reproduce plots provided in our papers and serve as code examples

- **exp_adoptercats.py** Neutral change among Roberts-style adopter categories. Vary the size of each category and the number of backward edges. contra Blythe &amp; Croft 2012

- **exp_corridory.py** Modeling the characteristic two-peak spread of NCS in the St. Louis Corridor and testing Friedman 2014's hypothesis

- **exp_singleg.py** Multi-clustered network structure yields S-curve change in populations of single-grammared speakers cf Yang 2009 cot-caught merger

- **exp_smallpops.py** Demonstrating the combined effect of small networks and single-grammar speakers on neutral (cf Kauhanen 2016) and advantaged change

## Model Files

These implement the model

- **applyT.py** Apply T learner transition matrix

- **calcP.py** Calculate P learner input matrix matrix

- **calcTfromextensions.py** Automatically create T matrix from grammar extensions

- **generateA.py** Generate A adjacency matrix

- **updateB.py** Update B grammar expression matrix

## Other Files

- **POC.py** Model proof of concept implemented in python/numpy

- **NBRecreate.py** Reproduces some results from Niyogi & Berwich 1997

- **V2_3param.py** Inputs for calcTfromextension.py** to reproduce Gibson & Wexler 1994