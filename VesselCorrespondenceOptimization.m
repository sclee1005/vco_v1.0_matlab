% VesselCorrespondenceOptimization.m : code of VCO (Vessel Correspondence Optimization) method 
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
% [2] Frangi filter (by D. Kroon) : 
%     https://kr.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
% [3] Frangi filter (code by D. Kroon, modified by Bingxiong Lin) : 
%   (slight improvement compared to D.Kroon frangi filter)
%     http://rpal.cse.usf.edu/project1/Vessel_Feature_Detection_Codes.zip
% [4] Fast marching method (by Gabriel Peyre) : 
%     https://kr.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching
% [5] MRF (trw-s, by Vladimir Kolmogorov) : 
%     https://github.com/aosokin/mrfMinimizerMex_trws_lbp
% [6] VLFeat : http://www.vlfeat.org/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bimg_tp1,all_v,bimg_tp1_post_processed] = ...
    VesselCorrespondenceOptimization(img_t,img_tp1,...
    bimg_t,p,pfrangi,pfm,psift)
%%% main function for VCO %%%
% inputs are pair of frames at t and t+1,(precomputed) vessel centerline of frame t, and
% parameters;
% output is the vessel centerline of t+1 frame

%%% summary of VCO subprocesses %%%
% 1. Global chafer matching
% 2. Point to point matching (get correspondence candidates for 'feature points' and 'vessel points')
% 3. MRF optimization (determine optimal correspondence candidates)
% 4. connect 'vessel points' to reform dense vessel centerline

%%% input parameter details %%%
% img_t : t frame vessel image
% img_tp1 : t+1 frame vessel image
% bimg_t : t frame mask image
% p : general parameters
% pfrangi : frangi vesselness parameters
% pfm : fast marching parameters
% psift : SIFT parameters
%%% output parameter details %%%
% bimg_tp1 : computed result of VCO befor post-processing
% all_v : corresponded sample points of t+1 frame
% bimg_tp1_post_processed : computed result of VCO after post-processing

%%% if there is no parameters, set up defaults paramters.
%%% mostly related to hierarchical matching %%%
if nargin < 4
    %%% general parameters
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
end

if nargin < 5
   %%% for frangi filter
    pfrangi = struct('FrangiScaleRange', [2 7], 'FrangiScaleRatio', 1,...
    'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 15, 'verbose',false,'BlackWhite',true);
end

if nargin < 6
    %%% for fast marching
    pfm.nb_iter_max = Inf; % parameters for fast marching
end

if nargin < 7
    %%% for SIFT
    psift.scale = 3.0; % scale factor for SIFT
    psift.magnif = 3; % magnitude
    psift.binSize = psift.magnif*psift.scale; % bin size
    psift.nnorm = 2; % l-2 norm
end
%%% if there is no parameters, set up defaults paramters.

% get to a image size
[nY,nX] = size(img_t);

% make frangi vesselness of tp1 frame
%%% either code [2] or [3] can be used, comment/uncomment as preferred
ivessel_tp1 = FrangiFilter2D(img_tp1*255, pfrangi); % use code of [2]
% ivessel_tp1 = FrangiFilter2DModifiedBingxiongLin(img_tp1, pfrangi); % use code of [3]

% bimg (binary mask of vessel centerline) of 't' frame
% preproceesing the binary mask image in case the input mask is not thinned mask
bimg_t = bwmorph(bimg_t,'thin',Inf);

% bimg of 't+1' frame
% construct estimate based on frangi vesselness for initial chamfer matching
bw = im2bw(ivessel_tp1,p.thre_ivessel);
bimg_tp1 = bwmorph(bw,'thin',Inf);

%% global chamfer matching
% for initial global translational motion estimation between vessels at frame t and t+1
% inputs are vessel centerlines at frame t and t+1 (here, coarse estimate of t+1 frame centerlines
% are computed from vesselness thresholding). output is translation vector (t_x, t_y) and
% translated t frame mask

% inputs
% bimg_t : t frame mask
% bimg_tp1 : t+1 frame mask which is frangi based mask
% p : parameters

% outputs
% gt_bimg_t : translated mask from t frame mask
% t_x, t_y : translated distance
[gt_bimg_t,t_x,t_y]=GlobalChaferMatching_run(bimg_t,bimg_tp1, p);

%% point-to-point matching(making putative matches)
% to generate corresponding point candidates in frame t+1 for points from
% frame t.
% Point2PointMatching funtion search correspondece point candidates using
% dense SIFT between t frame and t+1 frame

% inputs
% img_t : t frmae image
% gt_bimg_t : t frame mask
% img_tp1 : t+1 frame image
% ivessel_tp1 : t+1 frame vessleness
% t_x,t_y : translated distance from computed global chamfermatch
% p : paramters
% psift : sift parameters
% gc_canvas_img : t+1 frame image for visualiation
% save_path : save path
% fidx_tp1 : t+1 frame idex number

% outputs
% E : edges(piexel points) of t frame mask
% J : juntion points of t frame
% cell_cands : candidates of each subsampled points
% cell_cands_dists : computed distances from candidates of each subsampled points
[E,J,cell_cands,cell_cands_dists] = ...
    Point2PointMatching(img_t,gt_bimg_t,img_tp1,ivessel_tp1,t_x,t_y,p, psift);

%% MRF regularization
% to determine optimal correspondences using MRF optimization (TRW-S)

% compute costs for MRF optimization
[unaryCost,pairwiseCost,mapMat,all_coors,all_cands,cell_coors_to_all_coors] = ...
    ComputeCosts2(E,cell_cands,cell_cands_dists,p);
mapMat = sparse(mapMat);

% get each corresponded labels using MRF optimization
[labels,~] = mrfMinimizeMexModifiedForVCO(unaryCost,pairwiseCost,mapMat);
 
%% draw nodes & connect edges (basic post-processing)
% connect all output vco points of the t+1_th frame
% get all "feature" (bifurcations + crossings + end) point indices 

% inputs
% ivessel_tp1 : vesselness probability of t+1 frame
% nX,nY : image size
% E : edges of t frame(edge is subsampling point of t frame)
% J : juntion points of t frame
% all_coors : all of coordinats
% cell_coors_to_all_coors : segment based coordinates
% labels : each edges label which is resulted as MRF
% all_cands : all of candidates point of each edges point
% p : parameters

% outputs
% all_vessel_pt : all of result points of t+1 which are put in 1-D array
% all_v : result points of t+1 which are divided as segments
[all_vessel_pt,all_v] = ...
    MakeConnectedCenterlineFromSubsampledPts(ivessel_tp1,nX,nY,E,J,all_coors,...
    cell_coors_to_all_coors,labels,all_cands,p);

% for visualization
% all vessel points draw to image, and some process are added for smoothing
all_vessel_pt = unique(all_vessel_pt);
bimg_tp1 = false(nY,nX);
bimg_tp1(all_vessel_pt) = true;
bimg_tp1 = bwmorph(bimg_tp1,'dilate');
bimg_tp1 = bwmorph(bimg_tp1,'thin',Inf);
bimg_tp1 = bwmorph(bimg_tp1,'fill');
bimg_tp1 = bwmorph(bimg_tp1,'thin',Inf);
all_vessel_pt = find(bimg_tp1);

%% (additional) post-processing
% in case of angiograms, vessels are revealed as contrast agent flows in, which makes the vessels grow. 
% since this is not considered in VCO, here, we try to estimate the newly visible vessel regions.

% input
% ivessel : vesselness
% lidx : linear indices for vessels
% thre : threshold for 'ivessel', default 0.05
%
% output
% new_bimg : post-processed t+1 frame binary mask for considering extened vessels 
[new_bimg,~,~] = GrowVesselUsingFastMarching(ivessel_tp1,all_vessel_pt,p.thre_ivessel);

% copy post-processed t+1 frame mask and write
bimg_tp1_post_processed = new_bimg;
end