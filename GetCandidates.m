% GetCandidates.m : get candidates of 'feature points'
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

function [arr_cands,arr_cands_dists]= ...
    GetCandidates(img_t,E,tran_vec,d_tp1,idx_img_tp1,ivessel_tp1,p,psift)
% get candidates of 'feature points'

% inputs
% img_t : image of time 't', [0,1]
% E : a segment for which we will find correspondences, len*2, {{y,x}}
% tran_vec : translation by chamfer matching, (y,x)
% d_tp1 : (SIFT) descriptors of 't+1' frame, 128*(# of descriptors)
% idx_img_tp1 : index map for 'd_tp1', 'idx_img_tp1(y,x) = a' means that
%               the SIFT descriptor centered on (y,x) is at d_tp1(:,a)
% ivessel_tp1 : vesselness image of time 't', [0,1]
% p : parameters
% psift : SIFT parameters 
%
% outputs
% arr_cands : linear indices of candidates, (# of points)*(p.n_cands)
% arr_cands_dists : distances of candidates, (# of points)*(p.n_cands)
%
% coded by syshin(160229)

[nY,nX] = size(img_t);

tCC = sub2ind([nY,nX],E(:,1),E(:,2));

% sample points on the current vessel segment to reduce
% computational burdens
len = length(tCC);
samp_idx = [1:p.sampling_period:len];
if samp_idx(end) ~= len
    samp_idx = [samp_idx,len];
end
tCC = tCC(samp_idx);

npt = length(tCC);
arr_cands = zeros(npt,p.n_cands);
arr_cands_dists = inf(npt,p.n_cands);

for idx_pt = 1 : npt % for each point in this segment

    cpt = tCC(idx_pt);
    [cpt_yy,cpt_xx] = ind2sub([nY,nX],cpt);
    old_cpt_yy = cpt_yy-tran_vec(1);
    old_cpt_xx = cpt_xx-tran_vec(2);
    if old_cpt_yy > nY | old_cpt_yy < 1 | old_cpt_xx > nX | old_cpt_xx < 1
        continue;
    end
    old_cpt = sub2ind([nY,nX],old_cpt_yy,old_cpt_xx);

    fc = [old_cpt_xx;old_cpt_yy;psift.scale;0];
    [~,d_t] = vl_sift(single(img_t*255),'frames',fc) ;

    cand_dist = false(nY,nX);
    cand_dist(max(1,cpt_yy-p.img_bndry_th_m):min(nY,cpt_yy+p.img_bndry_th_m),...
        max(1,cpt_xx-p.img_bndry_th_m):min(nX,cpt_xx+p.img_bndry_th_m)) = true;
    cand_ivessel = ivessel_tp1>=p.thre_ivessel;
    cand = cand_ivessel&cand_dist;
    cand_idx = find(cand);
    dist_img = inf(nY,nX);
    for k = 1 : length(cand_idx)
        if idx_img_tp1(cand_idx(k))
            dist_img(cand_idx(k)) = norm(double(d_t)-d_tp1(:,idx_img_tp1(cand_idx(k))),psift.nnorm);
        end
    end
    [~,idx] = sort(dist_img(:),'ascend');
    cand_idx = idx(1:p.n_cands);
    arr_cands_dists(idx_pt,:) = dist_img(cand_idx);
    arr_cands(idx_pt,1:length(cand_idx)) = cand_idx';
    
end

end