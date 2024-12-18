# Bayesian_VariableSelection_GLMM

This repo contains the source code files that are associated with the paper [Stochastic Search Variable Selection for Bayesian Generalized Linear Mixed Effect Models](https://arxiv.org/pdf/2412.01084).  Variable selection remains a difficult problem, especially for generalized linear mixed models (GLMMs). While some frequentist approaches to simultaneously select joint fixed and random effects exist, primarily through the use of penalization, existing approaches for Bayesian GLMMs exist only for special cases, like that of logistic regression. In this work, we apply the Stochastic Search Variable Selection (SSVS) approach for the joint selection of fixed and random effects proposed in Yang et al. (2020) for linear mixed models to Bayesian GLMMs. We show that while computational issues remain, SSVS serves as a feasible and effective approach to jointly select fixed and random effects. We demonstrate the effectiveness of the proposed methodology to both simulated and real data. Furthermore, we study the role hyperparameters play in the model selection.


# Directory Introduction

## Codes
The directory contains all the source code files to implement our method on various GLMMs.

## simulation
The directory contains all the source code files required to run the simulation study in the paper.

## Real_data
The directory contains all the source code files required to run the real data analysis in the paper.
