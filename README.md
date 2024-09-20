# **Stochastic Modeling of Biofilm Formation with Quorum Sensing**

This repository contains the MATLAB code developed for the paper "**Stochastic Modeling of Biofilm Formation with Bacterial Quorum Sensing**," presented at the IEEE International Conference on Communications (ICC) 2023. The code implements a stochastic simulation method to model biofilm formation, focusing on the quorum sensing (QS) mechanism. The results generated from this code correspond to Figures 4-8 in the paper.

## **Abstract**

Bacteria often live in complex structures known as biofilms, composed of bacterial colonies and extracellular polymeric substances (EPS). Biofilms are associated with various issues, including infections and antibiotic resistance. This work proposes a stochastic model for biofilm formation using bacterial quorum sensing. The model captures biological processes like bacterial reproduction, autoinducer and EPS production, and their diffusion. A modified explicit tau-leap algorithm is used for simulations, adapted for the two-state QS mechanism. The model is validated using experimental data from *Pseudomonas putida* IsoF bacteria for autoinducer and bacterial concentration. The results show significant changes in EPS concentration before and after QS activation.

## **Code Overview**

This repository includes a single MATLAB code file that reproduces the results shown in Figures 4-8 of the paper. To model and simulate the formation of bacterial biofilms based on quorum sensing, the code implements a modified version of the tau-leap algorithm from the paper:

- Y. Cao, D. T. Gillespie, and L. R. Petzold, “Efficient step size selection for the tau-leaping simulation method,” The Journal of Chemical Physics, vol. 124, no. 4, p. 044109, 2006.

The modifications include downscaling autoinducer-related reaction rates and implementing a two-state quorum sensing mechanism. Additionally, the deterministic version of the model is solved using the MATLAB SimBiology Toolbox, which is adapted for variable threshold scenarios and numerical differentiation.

## **How to Run the Code**

- Simply run the provided MATLAB code file to generate the results corresponding to Figures 4-8 of the paper. Make sure you have the required MATLAB SimBiology Toolbox installed if you want to run the deterministic model comparisons. You can also run the code faster if you just uncomment the parfor() line and comment the for() line in the Modified Tau-Leap Algorithm section of the code. For using parfor(), Parallel Computing Toolbox is required.

## **Citation Requirement**

If you use or build upon this code in your research, or if this code is used for any academic work, publication, or research, proper attribution and citation of the paper is **required**. Please cite the paper below in any related publications, presentations, or derived works.

**Gulec, F., & Eckford, A. W., "Stochastic modeling of biofilm formation with bacterial quorum sensing," in ICC 2023 - IEEE International Conference on Communications, pp. 4470–4475, IEEE, 2023. https://ieeexplore.ieee.org/abstract/document/10278566. Arxiv link: [![arXiv](https://img.shields.io/badge/arXiv-2212.06269-b31b1b.svg)](https://arxiv.org/abs/2212.06269)**
