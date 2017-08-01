% MakeConnectedCenterlineFromSubsampledPts.m : Make connected center line 
% from subsampled points which are each labels from results ofMRF optimization.
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

function [all_vessel_pt,all_v] = ...
    MakeConnectedCenterlineFromSubsampledPts(ivessel_tp1,nX,nY,...
    E,J,all_coors,cell_coors_to_all_coors,labels,all_cands,p)
% Make connected center line from subsampled points which are each labels 
%                from results ofMRF optimization.

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

% coded by syshin(151104)

nE = size(E,1);
opt.method = 'discrete'; % set to discrete fast marching method

% detect junctions labeled as 'dummy'
bJ_all_pts = false(size(all_coors,1),1);

for j = 1 : nE
    if ismember(E{j}(1,:),J,'rows')
        bJ_all_pts(cell_coors_to_all_coors{j}(1)) = true;
    end
    if ismember(E{j}(end,:),J,'rows')
        bJ_all_pts(cell_coors_to_all_coors{j}(end)) = true;
    end
end
bJ_idx = find(bJ_all_pts);
all_joining_seg = {};
num_all_joining_seg = 0;
for j = 1 : length(bJ_idx)
    if labels(bJ_idx(j)) == p.n_all_cands+1
        %%% find segments joining at this junction
        joining_seg = [];
        for k = 1 : nE
            idx = find(cell_coors_to_all_coors{k}==bJ_idx(j));
            if idx
                if idx == 1
                    meet_pt = find(labels(cell_coors_to_all_coors{k})~=p.n_all_cands+1,1,'first');
                else
                    meet_pt = find(labels(cell_coors_to_all_coors{k})~=p.n_all_cands+1,1,'last');
                end
                if isempty(meet_pt)
                    continue;
                end
                joining_seg = [joining_seg;[k,meet_pt]];
            end
        end
        num_all_joining_seg = num_all_joining_seg +1 ;
        all_joining_seg{num_all_joining_seg} = joining_seg;   
        all_joining_junctions(num_all_joining_seg) = bJ_idx(j);
    end
end
% detect junctions labeled as 'dummy'

newE = {};
all_v = [];
all_vessel_pt = [];
cor_line_pt = [];
for j = 1 : nE
    temp_v = [];
    temp = false(nY,nX);
    t_seg = cell_coors_to_all_coors{j};
    len_t_seg = length(t_seg);
    st_pt = 0;
    ed_pt = 0;
    cum_seg_path = [];
    is_first = true;
    for k = 1 : len_t_seg
        t_idx = t_seg(k);
        t_label = labels(t_idx);
        if t_label < p.n_all_cands+1
            %%%
            pt1_y = all_coors(t_idx,1); pt1_x = all_coors(t_idx,2);
            pt2 = all_cands(t_idx,t_label);
            [pt2_y,pt2_x] = ind2sub([nY,nX],pt2);
            [t_path_x,t_path_y] = bresenham(pt1_x,pt1_y,pt2_x,pt2_y);
            cor_line_pt = [cor_line_pt;[t_path_y,t_path_x]];
            %%%
            if st_pt == 0
                st_pt = all_cands(t_idx,t_label);
                temp_v = [temp_v;st_pt];
            else
                ed_pt = all_cands(t_idx,t_label);
                temp_v = [temp_v;ed_pt];
            end
        else
            ed_pt = 0;
        end
        if st_pt ~= 0 && ed_pt ~= 0
            [st_pt_y,st_pt_x] = ind2sub([nY,nX],st_pt);
            [ed_pt_y,ed_pt_x] = ind2sub([nY,nX],ed_pt);
            % straight line
            %[t_path_x,t_path_y] = bresenham(st_pt_x,st_pt_y,ed_pt_x,ed_pt_y);
            % geodesic path
            pfm.end_points = [ed_pt_y;ed_pt_x];
            [D,~] = perform_fast_marchingModifiedForVCO(ivessel_tp1,[st_pt_y;st_pt_x], pfm);
            
            geo_path = compute_geodesicModifiedForVCO(D,[ed_pt_y;ed_pt_x],opt);
%             geo_path = compute_descret_geodes(D,[ed_pt_y;ed_pt_x],[]);
            geo_path = round(geo_path);
            [~, m, ~] = unique(geo_path','rows','first');
            geo_path = geo_path(:,sort(m))';
            geo_path = flipud(fliplr(geo_path));
            t_path_x = geo_path(:,1);
            t_path_y = geo_path(:,2);

            if is_first
                cum_seg_path = [cum_seg_path;[t_path_y,t_path_x]];
                is_first = false;
            else
                cum_seg_path = [cum_seg_path;[t_path_y(2:end),t_path_x(2:end)]];
            end
            st_pt = ed_pt;
            ed_pt = 0;
        end
    end
    newE{j} = cum_seg_path;
    if ~isempty(cum_seg_path)
        all_v = [all_v;temp_v];
        lidx = sub2ind([nY,nX],cum_seg_path(:,1),cum_seg_path(:,2));
        all_vessel_pt = [all_vessel_pt;lidx];
        temp(lidx) = true;
    end
end

% drawing for junctions labeled as 'dummy'
for j = 1 : num_all_joining_seg
    joining_seg = all_joining_seg{j};
    n_joining_seg = size(joining_seg,1);
    cum_path = [];
    for k = 1 : n_joining_seg          
        t_idx = cell_coors_to_all_coors{joining_seg(k,1)}(joining_seg(k,2));
        t_label = labels(t_idx);
        t_coor1 = all_cands(t_idx,t_label);
        [st_pt_y,st_pt_x] = ind2sub([nY,nX],t_coor1);
        for m = k+1 : n_joining_seg
            t_idx = cell_coors_to_all_coors{joining_seg(m,1)}(joining_seg(m,2));
            t_label = labels(t_idx);
            t_coor2 = all_cands(t_idx,t_label);
            [ed_pt_y,ed_pt_x] = ind2sub([nY,nX],t_coor2);

            % straight line
            %[t_path_x,t_path_y] = bresenham(st_pt_x,st_pt_y,ed_pt_x,ed_pt_y);
            % geodesic path
            pfm.end_points = [ed_pt_y;ed_pt_x];
            [D,~] = perform_fast_marchingModifiedForVCO(ivessel_tp1,[st_pt_y;st_pt_x], pfm);
            geo_path = compute_geodesicModifiedForVCO(D,[ed_pt_y;ed_pt_x],opt);
            geo_path = round(geo_path);
            [~, m, ~] = unique(geo_path','rows','first');
            geo_path = geo_path(:,sort(m))';
            geo_path = flipud(fliplr(geo_path));
            t_path_x = geo_path(:,1);
            t_path_y = geo_path(:,2);

            cum_path = [cum_path;[t_path_y(2:end-1),t_path_x(2:end-1)]];
        end
    end
    if isempty(cum_path)
        continue;
    end
    lidx = sub2ind([nY,nX],cum_path(:,1),cum_path(:,2));
    all_vessel_pt = [all_vessel_pt;lidx];    
end
% drawing for junctions labeled as 'dummy'

end

