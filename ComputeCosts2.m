% ComputeCosts2.m : compute costs for MRF optimization like unary turms, pair wise costs
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

function [unaryCost,pairwiseCost,mapMat,all_coors,...
    all_cands,cell_coors_to_all_coors] = ...
    ComputeCosts2(cell_coors,cell_cands,cell_cands_dists,p)
% compute costs for MRF optimization like unary turms, pair wise costs
% input
%
% cell_coors : cell for paths of each segment, nSegm*1 cell, nSegm is the number of segments
%           each segment has (nPT*d) values
% cell_cands : this contains candidate points per each (sampled) point
%           in each segment, nSegm*1 cell, nSegm is the number of segments
% cell_cands_dists : corresponding (unary) costs, nSegm*1 cell, 
%                   nSegm is the number of segments
% p : parameters
%
% output
%
% unaryCost : unary costs, nLabel*nNode
% pairwiseCost : pairwise costs, nEdge*1, each element is size of
%           (nLabel*nLabel)
% mapMat : mapping indices for 'pairwiseCost', 'mapMat(i.j) = k' means that
%       there is an edge between node i & j and its costs are in pairwiseCost{k}
% all_coors :
% all_cands : 
% cell_coors_to_all_coors : 
%
% coded by syshin (160130)

%% params
unary_thre = 400;
unary_trun_thre = 350;
dummy_unary_cost = unary_thre;
dummy_pairwise_cost1 = 75;
dummy_pairwise_cost2 = 4*dummy_pairwise_cost1;
dist_sigma = p.sampling_period/3;
alpha = 0.5;
beta = 0.5;

%% add nodes
nSegm = length(cell_coors);
nCand = size(cell_cands{1},2);
nLabel = nCand+1;
all_coors = [];
all_cands = [];
all_cands_dists = [];
all_coors = [all_coors;cell_coors{1}];
all_cands = [all_cands;cell_cands{1}];
all_cands_dists = [all_cands_dists;cell_cands_dists{1}];
cell_coors_to_all_coors = cell(nSegm,1);
cell_coors_to_all_coors{1} = [1:size(cell_coors{1},1)]';
for i = 2 : nSegm %%% for each segment
    temp = cell_coors{i};
    [mem,~] = ismember(temp,all_coors,'rows');
    all_coors = [all_coors;temp(~mem,:)];
    new_line = zeros(nnz(~mem),size(all_cands,2));
    new_line(:,1:nCand) = cell_cands{i}(~mem,:);
    all_cands = [all_cands;new_line];
    new_line = inf(nnz(~mem),size(all_cands_dists,2));
    new_line(:,1:nCand) = cell_cands_dists{i}(~mem,:);
    all_cands_dists = [all_cands_dists;new_line];

    [~,memIdx] = ismember(temp,all_coors,'rows');
    cell_coors_to_all_coors{i} = memIdx;
end
temp = all_cands_dists(:,nCand+1:size(all_cands_dists,2));
temp(temp==0) = inf;
all_cands_dists(:,nCand+1:size(all_cands_dists,2)) = temp;

nNode = size(all_coors,1);
% unary cost
% add a dummy label for nodes having no candidate
all_cands_dists(all_cands_dists>unary_thre) = inf;
all_cands_dists(all_cands_dists<=unary_thre&all_cands_dists>unary_trun_thre) = unary_trun_thre;
all_cands(all_cands_dists==inf) = 0;

% added for a redundancy check
temp_all_cands = all_cands;
temp_all_cands_dists = all_cands_dists;
all_cands = zeros(nNode,nCand);
all_cands_dists = inf(nNode,nCand);
for i = 1 : nNode
    [C,ia,~] = unique(temp_all_cands(i,:));
    if C(1) ~= 0
        unique_cands = temp_all_cands(i,ia);
        unique_cands_dists = temp_all_cands_dists(i,ia);
    else
        unique_cands = temp_all_cands(i,ia(2:end));
        unique_cands_dists = temp_all_cands_dists(i,ia(2:end));
    end
    all_cands(i,1:length(unique_cands)) = unique_cands;
    all_cands_dists(i,1:length(unique_cands_dists)) = unique_cands_dists;
end
% added for a redundancy check

unaryCost = [all_cands_dists';dummy_unary_cost*ones(1,nNode)];

%% add edges 
nEdge = 0;
nIntraEdge = 0;
pairwiseCost = {};
mapMat = zeros(nNode,nNode);
for i = 1 : nSegm %%% for each segment
    t_coors_to_all_coors = cell_coors_to_all_coors{i};
    for j = 1 : length(t_coors_to_all_coors)-1
        nEdge = nEdge+1;
        nIntraEdge = nIntraEdge+1;
        t_pCost1 = GetTruncatedPairwiseCost(all_coors(t_coors_to_all_coors(j),:),all_coors(t_coors_to_all_coors(j+1),:), ...
            all_cands(t_coors_to_all_coors(j),:),all_cands(t_coors_to_all_coors(j+1),:),dummy_pairwise_cost1,dummy_pairwise_cost2);
        
        t_pCost = t_pCost1; 
               
        pairwiseCost{nEdge} = t_pCost;
        mapMat(t_coors_to_all_coors(j),t_coors_to_all_coors(j+1))= nEdge;
    end
end

end

function cost_mat = GetTruncatedPairwiseCost(coor1,coor2,cands1,cands2,dummy_pairwise_cost1,dummy_pairwise_cost2)

nY = 512; nX = 512;
nCand = length(cands1);
nLabel = nCand+1;
[cands1_yy,cands1_xx] = ind2sub([nY,nX],cands1);
[cands2_yy,cands2_xx] = ind2sub([nY,nX],cands2);
diff1_y = cands1_yy-coor1(1);
diff1_x = cands1_xx-coor1(2);
diff2_y = cands2_yy-coor2(1);
diff2_x = cands2_xx-coor2(2);
diff_y = repmat(diff2_y,nCand,1)-repmat(diff1_y',1,nCand);
diff_x = repmat(diff2_x,nCand,1)-repmat(diff1_x',1,nCand);
diff = sqrt(diff_x.^2+diff_y.^2);

cost_mat = dummy_pairwise_cost1*ones(nLabel,nLabel);
diff = diff*10;
diff(diff>dummy_pairwise_cost2) = dummy_pairwise_cost2;
diff(cands1==0,:) = inf; % needless
diff(:,cands2==0) = inf; % needless
cost_mat(1:nCand,1:nCand) = diff;

end

function cost_mat = GetIntervalCost(pre_dist,dist_sigma,cands1,cands2,dummy_pairwise_cost1,dummy_pairwise_cost2)

nY = 512; nX = 512;
nCand = length(cands1);
nLabel = nCand+1;
[cands1_yy,cands1_xx] = ind2sub([nY,nX],cands1);
[cands2_yy,cands2_xx] = ind2sub([nY,nX],cands2);
diff_y = repmat(cands2_yy,nCand,1)-repmat(cands1_yy',1,nCand);
diff_x = repmat(cands2_xx,nCand,1)-repmat(cands1_xx',1,nCand);
diff = sqrt(diff_x.^2+diff_y.^2);

cost_mat = dummy_pairwise_cost1*ones(nLabel,nLabel);
fun = @(x) 1-exp(-(x-pre_dist).^2*0.5*dist_sigma^-2);
diff = fun(diff)*dummy_pairwise_cost2;

diff(cands1==0,:) = inf; % needless
diff(:,cands2==0) = inf; % needless
cost_mat(1:nCand,1:nCand) = diff;

end