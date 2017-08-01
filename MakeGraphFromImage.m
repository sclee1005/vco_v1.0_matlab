% MakeGraphFromImage.m : Get the graph information from a binary contour image
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

function [V,J,bJ,E,conn,mapMat] = MakeGraphFromImage(bimg)
% Get the graph information from a binary contour image
% input
% bimg : binary image of a contour
% output
% V : vertices, (numV*d) array
% J : juntion vertices
% bJ : junction indicator, 1 if the node is junction
% E : all edge paths (point array) constituting the graph,
%             cell(numE,1), numE is the number of edges
% conn : connectivity matrix, (numV*numV) array
% mapMat : mapping matrix for E, (numV*numV) array, mapMat(i,j) = 1 -> edge
%           path between node i & j is contained in E{1}
% coded by syshin(160129)

% parameter
patch_half_size = 5;

% get size
[nY,nX] = size(bimg);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

% junction & end point detection
branch_img = bwmorph(bimg,'branchpoints');
end_img = bwmorph(bimg,'endpoints');
[end_y,end_x] = find(end_img);

% added
dilated_branch_img = bwmorph(branch_img,'dilate');
CC = bwconncomp(dilated_branch_img,4);
numCC = CC.NumObjects;
branch_y = zeros(numCC,1);
branch_x = zeros(numCC,1);
branch_idx_img = zeros(nY,nX);
for j = 1 : numCC
    tCC = CC.PixelIdxList{j};
    branch_idx_img(tCC) = j;
    [yy,xx] = ind2sub([nY,nX],tCC);
    cy = round(mean(yy));
    cx = round(mean(xx));
    patch = ExtractPatchWithZeroPadding(bimg,[cy,cx],patch_half_size*2+1);
    [~,CPT] = bwdist(patch);
    [yy,xx] = ind2sub([patch_half_size*2+1,patch_half_size*2+1],CPT(patch_half_size+1,patch_half_size+1));
    cc = [yy,xx] + [cy-patch_half_size-1,cx-patch_half_size-1];
    branch_y(j) = cc(1);
    branch_x(j) = cc(2);
end
idx = find(ismember([end_y,end_x],[branch_y,branch_x],'rows'));
end_y(idx) = [];
end_x(idx) = [];
% added

J = [branch_y,branch_x];
V = [[branch_y,branch_x];[end_y,end_x]];
numV = size(V,1);
bJ = false(numV,1); bJ(1:size(J,1)) = true;
conn = false(numV,numV);
mapMat = zeros(numV,numV);

simg = bimg;
all_CC_pts = []; % junction CC
for j = 1 : numCC
    simg(CC.PixelIdxList{j}) = false;
    temp = false(nY,nX);
    temp(CC.PixelIdxList{j}) = true;
    idx = find(bimg&temp);
    all_CC_pts = [all_CC_pts;idx];
end
CC = bwconncomp(simg,8);
numCC = CC.NumObjects;
numE = numCC;
E = cell(numE,1);

for j = 1 : numCC
    tCC = CC.PixelIdxList{j};
    timg = uint8(bimg);
    timg(all_CC_pts) = 2;
    tE = [];
    idx = find(end_img(tCC));
    stV = 0; edV = 0;
    if length(idx) == 2 % case : end-end
        stIDX = tCC(1);
        timg(stIDX) = 0;
        [stY,stX] = ind2sub([nY,nX],stIDX);
        tE = [tE;[stY,stX]];
        bForwardFirst = true;
        bBackwardFirst = true;
        % forward path
        while true
            if bForwardFirst
                curY = stY; curX = stX;
                bForwardFirst = false;
            end
            temp = ExtractPatchWithZeroPadding(timg,[curY,curX],3);
            [incY,incX] = find(temp);
            nextY = curY+incY-(size(temp,1)-1); nextX = curX+incX-(size(temp,2)-1);
            if isempty(nextY)
                break;
            end
            curY = nextY(1); curX = nextX(1);
            timg(curY,curX) = 0;
            tE = [tE;[curY,curX]];
        end
        % backward path
        while true
            if bBackwardFirst
                curY = stY; curX = stX;
                bBackwardFirst = false;
            end
            temp = ExtractPatchWithZeroPadding(timg,[curY,curX],3);
            [incY,incX] = find(temp);
            nextY = curY+incY-(size(temp,1)-1); nextX = curX+incX-(size(temp,2)-1);
            if isempty(nextY)
                break;
            end
            curY = nextY; curX = nextX;
            timg(curY,curX) = 0;
            tE = [[curY,curX];tE];
        end
        stV = find(ismember(V,tE(1,:),'rows'));
        edV = find(ismember(V,tE(end,:),'rows'));
    elseif length(idx) == 1 % case : junction-end
        curIDX = tCC(idx);
        timg(curIDX) = 0;
        [curY,curX] = ind2sub([nY,nX],curIDX);
        edV = find(ismember(V,[curY,curX],'rows')); 
        tE = [tE;[curY,curX]];
        
        while true
            temp = ExtractPatchWithZeroPadding(timg,[curY,curX],3);
            [incY,incX] = find(temp);
            nextY = curY+incY-(size(temp,1)-1); nextX = curX+incX-(size(temp,2)-1);
            ii = find(timg(sub2ind([nY,nX],nextY,nextX))==2);
            if ~isempty(ii)
                if length(ii) > 1
                    dists = sum(abs([nextY(ii),nextX(ii)]-repmat([curY,curX],length(ii),1)),2);
                    [~,min_idx] = min(dists);
                    ii = ii(min_idx);
                end
                curY = nextY(ii); curX = nextX(ii);           
                branch_idx = branch_idx_img(curY,curX);  
                [t_path_x,t_path_y] = bresenham(curX,curY,branch_x(branch_idx),branch_y(branch_idx));
                stV = branch_idx;
                timg(sub2ind([nY,nX],t_path_y,t_path_x)) = 0;
                tE = [flipud([t_path_y,t_path_x]);tE];              
                break;
            end
            curY = nextY; curX = nextX;
            timg(curY,curX) = 0;
            tE = [[curY,curX];tE];
        end
    else % case : junction-junction
        stIDX = tCC(1);
        timg(stIDX) = 0;
        [stY,stX] = ind2sub([nY,nX],stIDX);
        tE = [tE;[stY,stX]];
        bForwardFirst = true;
        bBackwardFirst = true;
        
        % forward path
        while true
            if bForwardFirst
                curY = stY; curX = stX;
                bForwardFirst = false;
            end
            temp = ExtractPatchWithZeroPadding(timg,[curY,curX],3);
            [incY,incX] = find(temp);
            nextY = curY+incY-(size(temp,1)-1); nextX = curX+incX-(size(temp,2)-1);
            ii = find(timg(sub2ind([nY,nX],nextY,nextX))==2);
            if ~isempty(ii)
                if length(ii) > 1
                    dists = sum(abs([nextY(ii),nextX(ii)]-repmat([curY,curX],length(ii),1)),2);
                    [~,min_idx] = min(dists);
                    ii = ii(min_idx);
                end
                curY = nextY(ii); curX = nextX(ii);
                branch_idx = branch_idx_img(curY,curX);  
                [t_path_x,t_path_y] = bresenham(curX,curY,branch_x(branch_idx),branch_y(branch_idx));
                timg(sub2ind([nY,nX],t_path_y,t_path_x)) = 0;
                tE = [tE;[t_path_y,t_path_x]];
                break;
            end
            curY = nextY(1); curX = nextX(1);
            timg(curY,curX) = 0;
            tE = [tE;[curY,curX]];
        end
        % backward path
        while true
            if bBackwardFirst
                curY = stY; curX = stX;
                bBackwardFirst = false;
            end
            temp = ExtractPatchWithZeroPadding(timg,[curY,curX],3);
            [incY,incX] = find(temp);
            nextY = curY+incY-(size(temp,1)-1); nextX = curX+incX-(size(temp,2)-1);
            ii = find(timg(sub2ind([nY,nX],nextY,nextX))==2);
            if ~isempty(ii)
                if length(ii) > 1
                    dists = sum(abs([nextY(ii),nextX(ii)]-repmat([curY,curX],length(ii),1)),2);
                    [~,min_idx] = min(dists);
                    ii = ii(min_idx);
                end
                curY = nextY(ii); curX = nextX(ii);  
                branch_idx = branch_idx_img(curY,curX);  
                if j==7
                    aaa=0;
                end
                [t_path_x,t_path_y] = bresenham(curX,curY,branch_x(branch_idx),branch_y(branch_idx));
                timg(sub2ind([nY,nX],t_path_y,t_path_x)) = 0;
                tE = [flipud([t_path_y,t_path_x]);tE]; 
                break;
            end
            curY = nextY; curX = nextX;
            timg(curY,curX) = 0;
            tE = [[curY,curX];tE];
        end
        stV = find(ismember(V,tE(1,:),'rows'));
        edV = find(ismember(V,tE(end,:),'rows'));
    end
    E{j} = tE;
    conn(stV,edV) = true; conn(edV,stV) = true;
    mapMat(stV,edV) = j;
end

end

function patch = ExtractPatchWithZeroPadding(img,patch_center,patch_size)

[nY,nX] = size(img);
half_patch_size = floor(patch_size/2);
padded_img = zeros(nY+2*half_patch_size,nX+2*half_patch_size);
padded_img(half_patch_size+1:end-half_patch_size,half_patch_size+1:end-half_patch_size) = img;
ex_center = patch_center + half_patch_size;
patch = padded_img(ex_center(1)-half_patch_size:ex_center(1)+half_patch_size, ...
            ex_center(2)-half_patch_size:ex_center(2)+half_patch_size);
end