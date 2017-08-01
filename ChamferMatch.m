% ChamferMatch.m : compute minimum cost as chamfer matching
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

function matchImg = ChamferMatch(dtImg,tempImg,search_center,search_window_size)
% compute minimum cost as chamfer matching
% input
% dtImg : distance-map
% tempImg : template-map
% search_center : should be equal to the position (x,y) of 'left-top' corner of
% 'tempImg'
% search_window_size : default = 21 pixels
% output
% matchImg : match score map
%
% coded by syshin(151104)

[r1,c1] = size(dtImg);
[r2,c2] = size(tempImg);
matchImg = ones(r1,c1)*flintmax;

if nargin < 3
    umax = c1-c2;
    vmax = r1-r2;
    uc = floor(c2/2);
    vc = floor(r2/2);
    u_sp = 1;
    v_sp = 1;
else
    if nargin < 4
        search_window_size = 21;
    end
    half_search_window_size = floor(search_window_size/2);
    
    umax = min(c1-c2+1,search_center(1)+half_search_window_size);
    vmax = min(r1-r2+1,search_center(2)+half_search_window_size);
    uc = floor(c2/2);
    vc = floor(r2/2);
    u_sp = max(1,search_center(1)-half_search_window_size);
    v_sp = max(1,search_center(2)-half_search_window_size);
end

[~,addr] = GetAddressTable(tempImg,r1);

for u = u_sp : umax
    for v = v_sp : vmax
        offset = (u-1)*r1+v-1;
        score = sum(dtImg(offset+addr));
        matchImg(vc+v,uc+u)= score;
    end
end

end

function [nn,addr] = GetAddressTable(tempImg,stride)

% [r2,c2] = size(tempImg);
[rr,cc] = find(tempImg);
nn = length(rr);
addr = (cc-1)*stride+rr;
    
end