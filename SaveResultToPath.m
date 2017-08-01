% SaveResultToPath.m : For saving the VCO result
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

function [ ] = SaveResultToPath(img_tp1,bimg_tp1,all_v,...
    bimg_tp1_post_processed,save_path)
% For saving the VCO result
% inputs
% bimg_tp1 : result of VCO
% bimg_tp1_post_processed : post-processed result of VCO
% savePath : path for saving
% fidx : file index number

[nY, nX,~] = size(img_tp1);

%% post-processing before
% the original t+1 frame image copy to a canvas image
[y,x] = find(bimg_tp1);
idx = sub2ind([nY nX],y,x);
canvas_img = zeros(nY,nX,3);
canvas_img(:,:,1) = img_tp1;
canvas_img(:,:,2) = img_tp1;
canvas_img(:,:,3) = img_tp1;

canvas_img(idx) = false;
canvas_img(idx+nY*nX) = false;
canvas_img(idx+2*nY*nX) = true;
imwrite(canvas_img,strcat(save_path,sprintf('before_pp_rgb.bmp')));

canvas_img(all_v) = true;
canvas_img(all_v+nY*nX) = false;
canvas_img(all_v+2*nY*nX) = false;
imwrite(canvas_img,strcat(save_path,sprintf('before_pp_rgb_node.bmp')));


tp1_result_save_path = strcat(save_path,sprintf('before_pp_result.bmp'));
imwrite(bimg_tp1,tp1_result_save_path);

%% post-processing after

% t+1 frame correspondence sample points and a connected mask image draws 
% into a RGB-channel t+1 frame image for visualization and write them
[y,x] = find(bimg_tp1_post_processed);
idx = sub2ind([nY nX],y,x);

canvas_img(:,:,1) = img_tp1;
canvas_img(:,:,2) = img_tp1;
canvas_img(:,:,3) = img_tp1;

canvas_img(idx) = false;
canvas_img(idx+nY*nX) = false;
canvas_img(idx+2*nY*nX) = true;
imwrite(canvas_img,strcat(save_path,sprintf('after_pp_rgb.bmp')));

canvas_img(all_v) = true;
canvas_img(all_v+nY*nX) = false;
canvas_img(all_v+2*nY*nX) = false;
imwrite(canvas_img,strcat(save_path,sprintf('after_pp_rgb_node.bmp')));

tp1_pp_resut_save_path = strcat(save_path,sprintf('after_pp_result_pp.bmp'));
imwrite(bimg_tp1_post_processed,tp1_pp_resut_save_path);
end

