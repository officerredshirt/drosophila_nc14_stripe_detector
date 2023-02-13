# Semi-automated stripe detection pipeline for stage 5 fruit fly embryos

This MATLAB package implements a semi-automated processing pipeline for detecting and analyzing stripes of gene expression in early (stage 5) _Drosophila melanogaster_ embryos. An earlier version of this code was used in [Galupa et al. (2023)](https://doi.org/10.1016/j.devcel.2022.12.003).

The stripe detection algorithm was designed for the _eve_ gene but also functions for similarly distinct stripe patterns such as _ftz_ or _rho_. It relies on a segmented image with individually measured nuclear fluorescence levels.

Be forewarned that this is "research-grade" code that has not been optimized for general application. Please use functionality outside the main processing pipeline with caution.

<h2>Usage</h2>

We used this code on embryos fixed and processed with standard immunofluorescence stain or FISH protocols. The code for segmentation was adapted from the Garcia lab [mRNADynamics package](https://github.com/GarciaLab/mRNADynamics) and works on the DAPI channel. The code supports co-staining with built-in options for _sna_ or _fkh_ to normalize fluorescence and (in the case of _sna_) identify stripe width at a fixed point on the dorsoventral axis.

<h3>Directory structure</h3>
Files are assumed to be stored in [rootdir]/[proj]/[date] and to have the CZI file extension (if necessary this should be straightforward to modify since the code ultimately relies on the BioFormats library for image extraction). Results are saved in [rootdir]/processed_data/[proj]/[date] as MAT files with the same filename as the original image files.

<h3>Processing pipeline</h3>
The pipeline is invoked from the executable "main.m" and is outlined as follows:

1. Microscope images are loaded and Z-stacks are unpacked into individual TIFF files. As written it is assumed there is one CZI file per embryo, although the code can be modified to accommodate single CZI files with multiple embryos. The user-defined "histone channel" is used to generate a projectoin image that will be used for segmentation.
2. Embryos are processed sequentially. First the embryo is segmented from the histone channel projection. The user is then asked to manually curate the segmentation so as to identify outliers to exclude from analysis, or add/remove nuclei entirely.
3. Fluorescence levels are calculated for all identified nuclei as the average fluorescence intensity of pixels falling within each nucleus in a max projection of the fluroescent channel.
4. The embryo is autorotated to lie on an axis. The user manually checks the rotation and corrects if necessary. Note that answering "y" to the prompt to rotate 180 degrees will not allow the user to rotate the embryo again, in contrast to the "h" option to rotate 90 degrees, which will prompt the user for further rotations. In the following steps, stripes are assumed to be oriented vertically after embryo rotation.
5. The user is asked to draw a box around the region containing stripes. Usually this is the whole embryo, but if there are bright artifacts in the image it may be best to draw the box to exclude them. For best results, leave a gap between the edge of the detection box and the edge of each stripe. Once stripes are detected, the user is prompted to manually curate the identification as desired.
6. Manually identify the calibration (co-stain) region to which fluorescence levels in the stripe will be normalized. Skip this step if no co-stain or normalization is desired.
7. Stripe width is calculated and saved alongside other embryo data in a "fish_seg" object.

Once the data are generated for individual embryos, the user may optionally pool related embryos into a single "fish_pool" object that can be used for further analysis, for example, to calculate the mean stripe width across many samples. The user can then use the built-in plot functions for the "fish_pool" object or export the pooled data into a CSV file.
