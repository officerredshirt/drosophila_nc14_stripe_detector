# drosophila_nc14_stripe_detector
Semi-automated stripe detection pipeline for stage 5 fruit fly embryos
-----
This MATLAB package implements a semi-automated processing pipeline for detecting and analyzing stripes of gene expression in early (stage 5) _Drosophila melanogaster_ embryos. An earlier version of this code was used in [Galupa et al. 2023](https://doi.org/10.1016/j.devcel.2022.12.003). Be forewarned that this is "research-grade" code, so not all functions are implemented for optimal use in all cases.

Usage notes (to be updated further)
-----
We used this code on embryos fixed and processed with standard immunofluorescence stain or FISH protocols. The code for segmentation was adapted from the Garcia lab [mRNADynamics package](https://github.com/GarciaLab/mRNADynamics) and works on the DAPI channel.

The stripe detection algorithm was designed for the _eve_ gene but also functions for similarly distinct stripe patterns such as _ftz_ or _rho_. Note that stripes are assumed to be vertical after embryo rotation.
