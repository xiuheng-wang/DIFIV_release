%% README.txt -” A guide for the HSR package %%

Link to the datasets (Indian Pines and Salinas): http://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes
Link to the TensorLab toolbox for MATLAB: https://www.tensorlab.net
Link to the HySure package: https://github.com/alfaiate/HySure/tree/master/src

../demos contains:
	test_exp1.m: Plots figures of cost and SNR as a function of R1, R2, R3. You can choose to fix R2= R1 or R3. Data and figure saves automatically.

	test_exp2.m: Experiments on the reconstruction of spectral signature. You can plot all spectra, only portions, relative + normalized error. Uses data from text_exp1. Data and figure saves automatically.

	test_exp3.m: Experiments on metrics for different methods. You can pick tensor rank and multilinear ranks. Uses HySure results. Tables save automatically.

	test_exp4.m: Same experiments as test_exp1 but with HS-Pan fusion. Data and figure saves automatically.
RUN ONE OF THESE TO GET RESULTS 

../utils contains:
	crop.m: crops data to the specified dimensions

	gauss_kernel.m: Creates 1D gaussian filter with specified size and parameter sigma

	solveC_normal.m: used in tenRec to solve C with normal equation

	solveC_qr.m: used in tenRec to solve C with normal equation and QR factorization (for complexity purposes)

	spatial_deg.m: generates P1 and P2 with specified downsampling ratio

	spectral_deg.m: generates Pm from specified sensor (Landsat or Quickbird)

../methods contains:
	run_hosvd.m: runs HOSVD algorithm, returns estimated SRI and cell array of metrics

	run_sdf.m: runs HOSVD with optimization by the corresponding TensorLab model; returns estimated SRI and cell array of metrics

	stereo.m: STEREO algorithm as specified in the HSR paper (Sidiropoulos et. al.)

	tenRec.m: initialization for STEREO algorithm as specified in the HSR paper (Sidiropoulos et. al.)

../metrics contains:
	cc.m: cross-correlation

	ergas.m: global relative error

	nmse.m: normalized mean square error

	r_snr.m: signal-to-noise ration

	sam.m: spectral angle mapper
BETWEEN GROUNDTRUTH SRI AND ESTIMATE
