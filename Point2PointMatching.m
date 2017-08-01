% Point2PointMatching.m : to generate corresponding point candidates in frame t+1 for points from frame t.
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

function [E,J,cell_cands,cell_cands_dists] = ...
    Point2PointMatching(img_t,gt_bimg_t,img_tp1,ivessel_tp1,t_x,t_y,p,psift)
% to generate corresponding point candidates in frame t+1 for points from frame t.
% inputs
% img_t : t frame image
% gt_bimg_t : t frame binary center line mask
% img_tp1 : t+1 frame image
% ivessel_tp1 : t+1 frame vesselness map
% t_x,t_y : translated distance which is computed from global chafer matching
% p : parameters
% psift : parameters of sift

% ouputs
% E : edges(piexel points) of t frame mask
% J : juntion points of t frame
% cell_cands : candidates of each subsampled points
% cell_cands_dists : computed distances from candidates of each subsampled points
%% point-to-point matching(making putative matches)

% get image size
[nY,nX] = size(img_t);

% get informataion form t frame mask such as juntion points, edge points..
% input
% gt_bimg_t : t frame mask

% outpus
% J : junction points
% bJ : whether juntion point or no when all edge points
% E : edge points(= all of pixel points)
% mapMat : mapping information matrix
[~,J,bJ,E,~,mapMat] = MakeGraphFromImage(gt_bimg_t);

% DSIFT on 't+1' frame
smoothed_img_tp1 = vl_imsmooth(single(img_tp1*255), sqrt(psift.scale^2 -.25));
[f_tp1, d_tp1] = vl_dsift(smoothed_img_tp1, 'size', psift.binSize) ;

% convert to points format from 2D to 1D
d_tp1 = double(d_tp1);
f_tp1 = round(f_tp1);
idx_img_tp1 = zeros(nY,nX);
idx_img_tp1(sub2ind([nY,nX],f_tp1(2,:),f_tp1(1,:))) = [1:size(f_tp1,2)];

% junction matching 
nJ = size(J,1);
j_cands = zeros(nJ,p.n_junction_cands);
j_cands_dists = zeros(nJ,p.n_junction_cands);
num_j_flows = zeros(nJ,1);
j_flows = cell(nJ,1);

%%% for each junction
for j = 1 : nJ 
    % translation 
    cpt_yy = J(j,1); cpt_xx = J(j,2);
    old_cpt_yy = cpt_yy-t_y;
    old_cpt_xx = cpt_xx-t_x;

    % get SIFT desciptor
    fc = [old_cpt_xx;old_cpt_yy;psift.scale;0];
    [~,d_t] = vl_sift(single(img_t*255),'frames',fc) ;

    % find candidates use to SIFT descriptor
    cand_dist = false(nY,nX);
    cand_dist(max(1,cpt_yy-p.img_bndry_th_d):min(nY,cpt_yy+p.img_bndry_th_d),...
              max(1,cpt_xx-p.img_bndry_th_d):min(nX,cpt_xx+p.img_bndry_th_d)) = true;
    cand_ivessel = ivessel_tp1>=p.thre_ivessel;
    cand = cand_ivessel&cand_dist;
    cand_idx = find(cand);
    dist_img = inf(nY,nX);
    for k = 1 : length(cand_idx)
        if idx_img_tp1(cand_idx(k))
            dist_img(cand_idx(k)) = norm(double(d_t)-d_tp1(:,idx_img_tp1(cand_idx(k))),psift.nnorm);
        end
    end
    
    % sort distance and get candidates
    [~,idx] = sort(dist_img(:),'ascend');
    cand_idx = idx(1:p.n_junction_cands);
    j_cands(j,:) = cand_idx;
    j_cands_dists(j,:) = dist_img(cand_idx);

    %%% get displacement vectors for junctions using connected component
    %%% analysis & non-maximum supreesion
    tempBW = false(nY,nX);
    tempBW(cand_idx) = true;
    CC = bwconncomp(tempBW);
    nCC = CC.NumObjects;
    CCsize = zeros(nCC,1);
    Flows = [];
    for k = 1 : nCC
        tCC = CC.PixelIdxList{k};
        CCsize(k) = length(tCC);
        [~,idx] = min(dist_img(tCC));
        [endY,endX] = ind2sub([nY,nX],tCC(idx));
        tFlow = [endY,endX]-[cpt_yy,cpt_xx];
        Flows = [Flows;tFlow];
    end
    [~,idx] = sort(CCsize,'descend');
    n_sel = min(nCC,p.n_junction_cc_cands);
    sel_idx = idx(1:n_sel);
    num_j_flows(j) = n_sel;
    j_flows{j} = Flows(sel_idx,:);    
end

% local point matching 
nE = length(E);
cell_cands = cell(nE,1);
cell_cands_dists = cell(nE,1);

%%% for each segment
for j = 1 : nE 
    iTrial = 1;
    %%% original
    % get candidates all of current segment sample points
    [arr_cands,arr_cands_dists]= GetCandidates(img_t,E{j},[t_y,t_x],d_tp1,idx_img_tp1,ivessel_tp1,p,psift);
    
    % get candidates and distances
    npt = size(arr_cands,1);
    all_arr_cands = zeros(npt,p.n_all_cands);
    all_arr_cands_dists = inf(npt,p.n_all_cands);
    all_arr_cands(:,p.n_cands*(iTrial-1)+1:p.n_cands*iTrial) = arr_cands;
    all_arr_cands_dists(:,p.n_cands*(iTrial-1)+1:p.n_cands*iTrial) = arr_cands_dists;
    iTrial = iTrial + 1;
    
    %%% translated
    [stV,edV] = find(mapMat==j);
    is_stV_J = bJ(stV);
    is_edV_J = bJ(edV);
    if is_stV_J
        nFlows = num_j_flows(stV);
        for k = 1 : nFlows
            tFlow = j_flows{stV}(k,:);
            translated_E = E{j}+repmat(tFlow,size(E{j},1),1);
            if nnz(translated_E(:,1)<1|translated_E(:,2)<1|translated_E(:,1)>nY|translated_E(:,2)>nX)
                continue;
            end
            [arr_cands,arr_cands_dists]= GetCandidates(img_t,translated_E,[t_y+tFlow(1),t_x+tFlow(2)],d_tp1,idx_img_tp1,ivessel_tp1,p,psift);
            
            all_arr_cands(:,p.n_cands*(iTrial-1)+1:p.n_cands*iTrial) = arr_cands;
            all_arr_cands_dists(:,p.n_cands*(iTrial-1)+1:p.n_cands*iTrial) = arr_cands_dists;
            iTrial = iTrial + 1;
        end
    end
    if is_edV_J
        nFlows = num_j_flows(edV);
        for k = 1 : nFlows
            tFlow = j_flows{edV}(k,:);
            translated_E = E{j}+repmat(tFlow,size(E{j},1),1);
            if nnz(translated_E(:,1)<1|translated_E(:,2)<1|translated_E(:,1)>nY|translated_E(:,2)>nX)
                continue;
            end
            [arr_cands,arr_cands_dists]= GetCandidates(img_t,translated_E,[t_y+tFlow(1),t_x+tFlow(2)],d_tp1,idx_img_tp1,ivessel_tp1,p,psift);
            
            all_arr_cands(:,p.n_cands*(iTrial-1)+1:p.n_cands*iTrial) = arr_cands;
            all_arr_cands_dists(:,p.n_cands*(iTrial-1)+1:p.n_cands*iTrial) = arr_cands_dists;
            iTrial = iTrial + 1;
        end
    end
    cell_cands{j} = all_arr_cands;
    cell_cands_dists{j} = all_arr_cands_dists;
    
    ind = all_arr_cands(:);
    ind(ind==0) = [];
end

%%% for each segment
for j = 1 : nE 
    len = size(E{j},1);
    samp_idx = [1:p.sampling_period:len];
    if samp_idx(end) ~= len
        samp_idx = [samp_idx,len];
    end
    E{j} = E{j}(samp_idx,:);
end
end

