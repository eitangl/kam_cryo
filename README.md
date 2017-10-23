# Ab initio Kam's method for cryo-EM data

## Prerequisites:
1. [ASPIRE](http://spr.math.princeton.edu/)
2. [Manopt](http://manopt.org/)
3. [Spherical Harmonic Transform toolbox by Archontis Politis](https://www.mathworks.com/matlabcentral/fileexchange/43856-real-complex-spherical-harmonic-transform--gaunt-coefficients-and-rotations)

## Installation:

1. Download and install ASPIRE, add to path by running the initpath.m script.
2. Download and install Manopt, add to path by running the importmanopt.m scriot.
3. Download and extract the Spherical Harmonics toolbox, add its folder to the path.
4. Download and extract this package [Warning: It includes a ~1GB file.]
5. (Optional) Run gen_proj_ASPIRE on a given volume to generate a number of uniformly randomly distrinuted projections to test synthetic datasets.

## Running a test example:

We provide the file 'empiar_10107_200K.mat' including the denoised images and estiamted covariance matrix from the EMPIAR-10107 dataset. 
You may run 'raw_data_rec_script.m' to process this file and obtain an ab initio model.

## Main functions:

1. 'raw_data_process_script.m': Processes raw experimental data from a star file, estimates covariance matrix and denoises a given number of images.
2. 'raw_data_rec_script.m': Produces *ab inito* reconstructions from the output of 'raw_data_process_script.m'.
3. 'run_exp_cryo.m': Gets a clean datasets, generates synthetic data with specified parameters and produces a reconstruction.

In case of issues or questions, please email Eitan (eitanl@math.princeton.edu)
