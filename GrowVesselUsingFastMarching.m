% GrowVesselUsingFastMarching.m : post-processing to estimate the newly visible vessel regions
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

function [new_bimg,new_lidx,app_lidx] = ...
    GrowVesselUsingFastMarching(ivessel,lidx,thre)
% after VCO, it was not considered for growing vessels
% so we get grown vessels, make this post-processing function

% input
% ivessel : vesselness
% lidx : linear indices for vessels
% thre : threshold for 'ivessel', default 0.05
%
% output
% new_bimg : binary mask for a new vessel
% new_lidx : linear indices for a new vessels
% app_lidx : linear indices of appened parts
%
% coded by syshin(160305)

if nargin < 3 
    thre = 0.05;
end
opt.method = 'discrete'; % set to discrete fast marching method

verbose=false;
IS3D=false;

[nY,nX] = size(ivessel);

% Convert double image to logical
Ibin=ivessel>=thre;
CC = bwconncomp(Ibin);
numCC = CC.NumObjects;
bROI = false(numCC,1);
Ibin = false(nY,nX);
for i = 1 : numCC
    Lia = ismember(CC.PixelIdxList{i},lidx);
    if ~isempty(find(Lia))
        bROI(i) = true;
        Ibin(CC.PixelIdxList{i}) = true;
    end
end

% Distance to vessel boundary
BoundaryDistance=getBoundaryDistance(Ibin,IS3D);


if(verbose)
    disp('Distance Map Constructed');
end

% Get maximum distance value, which is used as starting point of the
% first skeleton branch
[~,maxD]=maxDistancePoint(BoundaryDistance,Ibin,IS3D);

% Make a fastmarching speed image from the distance image
SpeedImage=(BoundaryDistance/maxD).^4;
SpeedImage(SpeedImage==0)=1e-10;

% Skeleton segments found by fastmarching
SkeletonSegments=cell(1,1000);

% Number of skeleton iterations
itt=0;

[yy,xx] = ind2sub([nY,nX],lidx);
SourcePoint = [yy';xx'];
 
while(true)
    if(verbose),
        disp(['Find Branches Iterations : ' num2str(itt)]);
    end

    % Do fast marching using the maximum distance value in the image
    % and the points describing all found branches are sourcepoints.
%     [T,Y] = msfm(SpeedImage, SourcePoint, false, false); %%origin code
    
    %%kjn
    [T,S] = perform_fast_marchingModifiedForVCO(ivessel,SourcePoint);
    
    tmp = zeros(512,512);
    for i=1:size(SourcePoint,2)
        tmp(SourcePoint(1,i),SourcePoint(2,i)) =1;
    end
    
    Y = double(bwdist(tmp));
    %%kjn

    % Trace a branch back to the used sourcepoints
    StartPoint=maxDistancePoint(Y,Ibin,IS3D);
    
    
    %ShortestLine=shortestpath(T,StartPoint,SourcePoint,1,'rk4'); %%origin code
    
    %%kjn
    ShortestLine = round(compute_geodesicModifiedForVCO(T,StartPoint,opt)); %% kjn
    ShortestLine = ShortestLine';

    ShortestLine = unique(ShortestLine,'rows','first'); %% kjn
    %%kjn
    
    
    % Calculate the length of the new skeleton segment
    linelength=GetLineLength(ShortestLine,IS3D);
        
    % Stop finding branches, if the lenght of the new branch is smaller
    %%
    % 
    %   for x = 1:10
    %       disp(x)
    %   end
    % 
    % then the diameter of the largest vessel
    if(linelength<maxD*2), break; end;
    
    % Store the found branch skeleton
    itt=itt+1;
    SkeletonSegments{itt}=ShortestLine;
    
    % Add found branche to the list of fastmarching SourcePoints
    SourcePoint=[SourcePoint ShortestLine'];
end
SkeletonSegments(itt+1:end)=[];

lidx_app = [];
if ~isempty(SkeletonSegments)

S=OrganizeSkeleton(SkeletonSegments,IS3D);
if(verbose),
    disp(['Skeleton Branches Found : ' num2str(length(S))]);
end


for i=1:length(S)
    L=round(S{i});
    [b, m, n] = unique(L','rows','first');
    L = L(:,sort(m));
    S{i} = L;
    lidx_app = [lidx_app;sub2ind([nY,nX],L(:,1),L(:,2))];
end

end

% % Display the skeleton
% figure, imshow(Ibin); hold on;
% for i=1:length(S)
% L=S{i};
% plot(L(:,2),L(:,1),'-','Color',rand(1,3));
% end

app_lidx = lidx_app;
new_lidx = [lidx;lidx_app];
new_bimg = false(nY,nX);
new_bimg(new_lidx) = true;

end

function ll=GetLineLength(L,IS3D)
if(IS3D)
    dist=sqrt((L(2:end,1)-L(1:end-1,1)).^2+ ...
              (L(2:end,2)-L(1:end-1,2)).^2+ ...
              (L(2:end,3)-L(1:end-1,3)).^2);
else
    dist=sqrt((L(2:end,1)-L(1:end-1,1)).^2+ ...
              (L(2:end,2)-L(1:end-1,2)).^2);
end
ll=sum(dist);
end

    
function S=OrganizeSkeleton(SkeletonSegments,IS3D)
n=length(SkeletonSegments);
if(IS3D)
    Endpoints=zeros(n*2,3);
else
    Endpoints=zeros(n*2,2);
end
l=1;
for w=1:n
    ss=SkeletonSegments{w};
    l=max(l,length(ss));
    Endpoints(w*2-1,:)=ss(1,:); 
    Endpoints(w*2,:)  =ss(end,:);
end
CutSkel=spalloc(size(Endpoints,1),l,10000);
ConnectDistance=2^2;

for w=1:n
    ss=SkeletonSegments{w};
    ex=repmat(Endpoints(:,1),1,size(ss,1));
    sx=repmat(ss(:,1)',size(Endpoints,1),1);
    ey=repmat(Endpoints(:,2),1,size(ss,1));
    sy=repmat(ss(:,2)',size(Endpoints,1),1);
    if(IS3D)
        ez=repmat(Endpoints(:,3),1,size(ss,1));
        sz=repmat(ss(:,3)',size(Endpoints,1),1);
    end
    if(IS3D)
        D=(ex-sx).^2+(ey-sy).^2+(ez-sz).^2;
    else
        D=(ex-sx).^2+(ey-sy).^2;
    end
    check=min(D,[],2)<ConnectDistance;
    check(w*2-1)=false; check(w*2)=false;
    if(any(check))
        j=find(check);
        for i=1:length(j)
            line=D(j(i),:);
            [foo,k]=min(line);
            if((k>2)&&(k<(length(line)-2))), CutSkel(w,k)=1; end
        end
    end
end

pp=0;
for w=1:n
    ss=SkeletonSegments{w};
    r=[1 find(CutSkel(w,:)) length(ss)];
    for i=1:length(r)-1
        pp=pp+1;
        S{pp}=ss(r(i):r(i+1),:);
    end
end
end

function BoundaryDistance=getBoundaryDistance(I,IS3D)
% Calculate Distance to vessel boundary

% Set all boundary pixels as fastmarching source-points (distance = 0)
if(IS3D),S=ones(3,3,3); else S=ones(3,3); end
B=xor(I,imdilate(I,S));
ind=find(B(:));
if(IS3D)
    [x,y,z]=ind2sub(size(B),ind);
    SourcePoint=[x(:) y(:) z(:)]';
else
    [x,y]=ind2sub(size(B),ind);
    SourcePoint=[x(:) y(:)]';
end

% Calculate Distance to boundarypixels for every voxel in the volume
SpeedImage=ones(size(I));
% BoundaryDistance = msfm(SpeedImage, SourcePoint, false, true); %% origin code

BoundaryDistance = bwdist(B,'euclidean');%%kjn

% Mask the result by the binary input image
BoundaryDistance(~I)=0;
end

function [posD,maxD]=maxDistancePoint(BoundaryDistance,I,IS3D)
% Mask the result by the binary input image
BoundaryDistance(~I)=0;

% Find the maximum distance voxel
[maxD,ind] = max(BoundaryDistance(:));
if(~isfinite(maxD))
    error('Skeleton:Maximum','Maximum from MSFM is infinite !');
end

if(IS3D)
    [x,y,z]=ind2sub(size(I),ind); posD=[x;y;z];
else
    [x,y]=ind2sub(size(I),ind); posD=[x;y];
end
end