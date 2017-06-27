# BSD
Blind Spare Deconvolution for inferring spike trains from fluorescence recordings


Author: Jérôme Tubiana (jertubiana@gmail.com)

Affiliations:  
- Laboratoire de Physique Théorique, Ecole Normale Supérieure, PSL, Paris.
-  Laboratoire Jean Perrin, UPMC, PSL, Paris.

Copyright: MIT License, 2017

Reference Article: http:/biorxiv etc.

## Installation

Requirements: Matlab with signal processing & optimisation toolboxes.
Installation: Add the folder to matlab path:

```
addpath(genpath(‘path_to_BSD_folder/BSD’))
```

## Contents
- BSD.m: the main function.
- pBSD.m: the parallel version for handling multiple fluorescence traces.
- BSD_theoretical_accuracy.m: A function for evaluating theoretical limits (precision-recall & temporal accuracy)
- BSD & pBSD main Subroutines:
    - BSD_initialization.m: Initial generative model parameters estimation form fluorescence trace.
    - BSD_deconvolution.m: Sparse deconvolution from fluorescence trace and generative model parameters
    - BSD_parameter_estimation.m: Refinement of generative model parameters from inferred spikes & fluorescence trace.
- utilities/… : other subroutines.
- examples/.. : folder containing examples of usage of the main functions.

## Usage
Refer to the scripts in the examples folder for usage instruction.
