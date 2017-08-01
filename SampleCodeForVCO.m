% SampleCodeForVCO.m : sample code to run VCO (Vessel Correspondence Optimization) method 
%
% VCO method described in
% [1] Seung Yeon Shin, Soochahn Lee, Kyoung Jin Noh, Il Dong Yun, Kyoung Mu Lee:
%     Extraction of Coronary Vessels in Fluoroscopic X-Ray Sequences Using 
%     Vessel Correspondence Optimization. 
%     MICCAI (3) 2016: 308-316
%
% code by : Seung Yeon Shin, Soochahn Lee, Kyoung Jin Noh
% Last updated : 2017-07-18
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Links for 3rd party code
% 1. Frangi filter - by D. Kroon
% https://kr.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
% 2. Frangi filter (code by D. Kroon, modified by Bingxiong Lin)
%   (slight improvement compared to D.Kroon frangi filter)
% http://rpal.cse.usf.edu/project1/Vessel_Feature_Detection_Codes.zip
% 3. Fast marching method - by Gabriel Peyre
% https://kr.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching
% 4. MRF (trw-s) - by Vladimir Kolmogorov
% https://github.com/aosokin/mrfMinimizerMex_trws_lbp
% 5. VLFeat
% http://www.vlfeat.org/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%% VCO parameters %%%%%%%%%%%%%%%%%%%%%
%%% mostly related to hierarchical matching %%%
p.use_gc = true; % use global chamfer matching or not
p.n_junction_cc_cands = 2; % N_k (in [1]) - # of candidates for 'feature points' 
p.n_cands = 5; % N_l (in [1]) - # of candidates for 'vessel points' for vessel segment
p.n_all_cands = (1+2*p.n_junction_cc_cands)*p.n_cands; % N_p (in [1]) - maximum # of total candidates
p.sampling_period = 5; % to control sampling of 'vessel points' on vessel segment (higher value = more sparse points)
%%% internal parameters %%%
p.n_junction_cands = 10; % # of 'feature points' candidates examined before non-max suppression
p.thre_ivessel = 0.05; % threshold for frangi vesselness used for t+1 frame 
p.img_bndry_th_d = 50; % image boundary threshold : to exclude points near image boundary (description)
p.img_bndry_th_m = 10; % image boundary threshold : to exclude points near image boundary (matching)
%%% for frangi filter
pfrangi = struct('FrangiScaleRange', [2 7], 'FrangiScaleRatio', 1,...
'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 15, 'verbose',false,'BlackWhite',true);
%%% for fast marching
pfm.nb_iter_max = Inf; % parameters for fast marching
%%% for SIFT
psift.scale = 3.0; % scale factor for SIFT
psift.magnif = 3; % magnitude
psift.binSize = psift.magnif*psift.scale; % bin size
psift.nnorm = 2; % l-2 norm
%%%%%%%%%%%%%%%%%%

%% data loading
% load frame t image file
t_img_path = 'test_img/img_t.bmp';
img_t = imread(t_img_path); % get a t frmae image 
img_t = im2double(img_t(:,:,1));
% load frame t mask file
t_mask_path = 'test_img/mask_t.bmp';
bimg_t = imread(t_mask_path); % get a t frame mask image
bimg_t = im2double(bimg_t(:,:,1));
% load frame t+1 image file
tp1_img_path = 'test_img/img_tp1.bmp';
img_tp1 = imread(tp1_img_path); % get a t+1 frame image
img_tp1 = im2double(img_tp1(:,:,1));

%% run VCO
% MAIN FUNCTION FOR VCO
% inputs
% t frame image, t+1 frame image, t frame mask, t+1 frame mask, parameters

% ouputs
% t+1 frame mask, t+1 frame mask 1-D array, post-processed t+1 frame mask
[bimg_tp1,all_v,bimg_tp1_post_processed] = VesselCorrespondenceOptimization(img_t,img_tp1,bimg_t,...
    p,pfrangi,pfm,psift);
%% save resut image write function
save_path = 'res/'; % result mask path for saving 
SaveResultToPath(img_tp1,bimg_tp1,all_v,bimg_tp1_post_processed,save_path)
