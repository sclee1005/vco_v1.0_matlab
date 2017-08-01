% GlobalChaferMatching_run.m : run Global chamfer matching
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

function [gt_bimg_t, t_x, t_y ] = GlobalChaferMatching_run(bimg_t,bimg_tp1, p )
%% global chamfer matching

% inputs 
% bimg_t : a t frame binary mask image
% bimg_tp1 : a t+1 frame binary mask image
% p : parameters
% save_path : save path
% fidx_tp1 : file index
% gc_canvas_img  : a t+1 frame RGB channel image for visualization

% ouputs
% gt_bimg_t : t frame mask is translated by global chamfermatching result
% t_x, t_y : tanslated distance

% get image size
[nY,nX] = size(bimg_t);

% do you compute global chamfer or not
if p.use_gc  
    % get distance transform
    dt_tp1 = bwdist(bimg_tp1);
    
    % croping t frame mask
    [rows,cols] = find(bimg_t);
    minX = min(cols);
    maxX = max(cols);
    minY = min(rows);
    maxY = max(rows);
    if rem(maxX-minX,2) == 1
        if maxX < nX
            maxX = maxX+1;
        else
            minX = minX-1;
        end
    end
    if rem(maxY-minY,2) == 1
        if maxY < nY
            maxY = maxY+1;
        else
            minY = minY-1;
        end
    end
    minY = max(1,minY);maxY = min(nY,maxY);
    minX = max(1,minX);maxX = min(nX,maxX);
   
    cropped_bimg_t = bimg_t(minY:maxY,minX:maxX);

    % do chamfer matching between cropped t frame mask and t+1 frame mask
    matchImg = ChamferMatch(dt_tp1,cropped_bimg_t);

    % bundary costs make to infinit value
    cpt_yy = round((minY+maxY)/2);
    cpt_xx = round((minX+maxX)/2);
    cand = false(nY,nX);
    cand(max(1,cpt_yy-p.img_bndry_th_d):min(nY,cpt_yy+p.img_bndry_th_d),...
         max(1,cpt_xx-p.img_bndry_th_d):min(nX,cpt_xx+p.img_bndry_th_d)) = true;
    idx = find(~cand);
    matchImg(idx) = flintmax;
    
    % find to a minimum cost position
    [minval,idx] = min(matchImg(:));
    [yy,xx] = ind2sub([nY,nX],idx);
    hy = floor(size(cropped_bimg_t,1)/2);
    hx = floor(size(cropped_bimg_t,2)/2);

    % compute translated distance
    t_x = xx-hx-minX;
    t_y = yy-hy-minY;
    if minval == flintmax
        t_x = 0;
        t_y = 0;
    end
    
    % make translated t frame mask and write
    [byy,bxx] = find(bimg_t);
    byy = byy+t_y;
    bxx = bxx+t_x;
    gt_bimg_t = false(nY,nX);
    gt_bimg_t(sub2ind([nY,nX],byy,bxx)) = true;
else
    % if you do not global chafer matching
    % traslated distance is 0 
    gt_bimg_t = bimg_t;
    t_x = 0; t_y = 0;
end


end

