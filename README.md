# vco_v1.0_matlab
Matlab implementation (ver. 1.0) of the vessel correspondence optimization (VCO) algorithm described in:

    Seung Yeon Shin, Soochahn Lee, Kyoung Jin Noh, Il Dong Yun, Kyoung Mu Lee:
    Extraction of Coronary Vessels in Fluoroscopic X-Ray Sequences Using 
    Vessel Correspondence Optimization. 
    MICCAI (3) 2016: 308-316

**************************************************************
Brief description of important files
- SampleCodeForVCO.m: sample code demonstrating how to call the VesselCorrespondenceOptimization() function
- VesselCorrespondenceOptimization.m: main file for the VCO algorithm
- GlobalChaferMatching_run.m: global chamfer matching from centerlines of frame t to frame t+1
- Point2PointMatching.m: to generate corresponding point candidates in frame t+1 for points from frame t
- ComputeCosts2: compute MRF unary and pairwise energy costs
- mrfMinimizeMex_syshin.m (built after compiling buildMrfMinimizeMexModifiedForVCO.m in /third_party/mrfMinimizerMex_trws_lbp-master): perform MRF energy minimization to compute the optimal point correspondences
- GrowVesselUsingFastMarching.m: post-processing for data such as angiograms where visible vessel regions may be growing

**************************************************************
How to run:
1. Extract files
2. Download, copy, and where needed compile the following 3rd party code to the corresponding folders in the /third_party subfolder. (Make sure to add /third_party and all subfolders to Matlab path)
3. Run SampleCodeForVCO to obtain VCO results for provided sample images. Modify corresponding file paths in SampleCodeForVCO.m to run VCO on different images.
**************************************************************
Links for 3rd party code
- [1] Frangi filter : https://kr.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
  Included within the folder is a modified implementation of the Frangi filter which is a slight improvement compared to the above implementation - link for original code :  http://rpal.cse.usf.edu/project1/Vessel_Feature_Detection_Codes.zip
- [2] Fast marching method* : https://kr.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching
  * to compile, instead of running the original compile_mex.m file, run compile_mexModifiedForVCO.m in order for our modifications to take effect.
- [3] MRF* : https://github.com/aosokin/mrfMinimizerMex_trws_lbp 
  * to compile, instead of running the original buildMrfMinimizeMex.m file, run buildMrfMinimizeMexModifiedForVCO.m in order for our modifications to take effect.
- [4] VLFeat : http://www.vlfeat.org/
**************************************************************
If you have any question or recommendation regarding this code, please contact:
Seung Yeon Shin (syshin@snu.ac.kr), 
Seoul National University, Seoul, Republic of Korea
Soochahn Lee (sclsch@sch.ac.kr),
Kyoung Jin Noh (yellowd91@gmail.com)
Soonchunhyang University, Asan, Republic of Korea


More information about this project can be found at:
http://cv.snu.ac.kr/research/vco/


Seung Yeon Shin, Soochahn Lee, Kyoung Jin Noh
July 31, 2017

