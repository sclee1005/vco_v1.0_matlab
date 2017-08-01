function [outIm, whatScale, Direction, outDoH] = FrangiFilter2DModifiedBingxiongLin(I, options)
% This function FRANGIFILTER2D uses the eigenvectors of the Hessian to
% compute the likeliness of an image region to vessels, according
% to the method described by Frangi:2001 (Chapter 2).
%
% [J,Scale,Direction] = FrangiFilter2D(I, Options)
%
% inputs,
%   I : The input image (vessel image)
%   Options : Struct with input options,
%       .FrangiScaleRange : The range of sigmas used, default [1 8]
%       .FrangiScaleRatio : Step size between sigmas, default 2
%       .FrangiBetaOne : Frangi correction constant, default 0.5
%       .FrangiBetaTwo : Frangi correction constant, default 15
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.
%       .verbose : Show debug information, default true
%
% outputs,
%   J : The vessel enhanced image (pixel is the maximum found in all scales)
%   Scale : Matrix with the scales on which the maximum intensity 
%           of every pixel is found
%   Direction : Matrix with directions (angles) of pixels (from minor eigenvector)   
%               atan2: [-pi pi]
%   DoH: determinant of Hessian: lambda1*lambda2
%
% Example,
%   I=double(imread ('vessel.png'));
%   Ivessel=FrangiFilter2D(I);
%   figure,
%   subplot(1,2,1), imshow(I,[]);
%   subplot(1,2,2), imshow(Ivessel,[0 0.25]);
%
% Written by Marc Schrijver, 2/11/2001
% Re-Written by D.Kroon University of Twente (May 2009)
% modified by Bingxiong Lin, Jan. 2014
flag_calc_original_frangi = false;
flag_use_original_frangi = false;
if max(max(I))<=1
    I=I*255; %[0, 255] required by frangifilter
end
% defaultoptions = struct('FrangiScaleRange', [2 7], 'FrangiScaleRatio', 2.5, 'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 15, 'verbose',true,'BlackWhite',true);
defaultoptions = struct('FrangiScaleRange', [2 7], 'FrangiScaleRatio', 2.5, 'FrangiBetaOne', 0.5, 'FrangiBetaTwo', 15, 'verbose',false,'BlackWhite',true);
%be reminded scale in vesselness is actually the sigma (standard deviation) of gaussian kernal
% Process inputs
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('FrangiFilter2D:unknownoption','unknown options found');
    end
end

sigmas=options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2);
sigmas = sort(sigmas, 'ascend');

beta  = 2*options.FrangiBetaOne^2;
c     = 2*options.FrangiBetaTwo^2;

% Make matrices to store all filterd images
ALLfiltered=zeros([size(I) length(sigmas)]);
ALLangles=zeros([size(I) length(sigmas)]);

% Frangi filter for all sigmas
for i = 1:length(sigmas),
    % Show progress
    if(options.verbose)
        disp(['Current Frangi Filter Sigma: ' num2str(sigmas(i)) ]);
    end
    
    % Make 2D hessian
    [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));
    
    % Correct for scale
    Dxx = (sigmas(i)^2)*Dxx;
    Dxy = (sigmas(i)^2)*Dxy;
    Dyy = (sigmas(i)^2)*Dyy;
   
    % Calculate (abs sorted) eigenvalues and vectors
    [Lambda1,Lambda2,Ix,Iy, mu1, mu2]=eig2imageModifiedBingxiongLin(Dxx,Dxy,Dyy);%abs(Lambda1)<abs(Lambda2), mu1<mu2    
    
    % Compute the direction of the minor eigenvector
    angles = atan2(Iy,Ix);
    % Compute some similarity measures
    Lambda2(Lambda2==0) = eps;
    Rb = (Lambda1./(Lambda2)).^2;
    if(options.BlackWhite)
        mu1(mu1<0) = 0;
        mu2(mu2<0) = 0;
        %mu1 = mu1.*mu2;
    else
        mu1(mu1>0) = 0;
        mu2(mu2>0) = 0;
    end
    S2 = mu1.^2 + mu2.^2;  %
    % Compute the output image
    Ifiltered =  (ones(size(I))-exp(-S2/c));
    ALLfiltered(:,:,i) = Ifiltered;
    ALLangles(:,:,i) = angles;
    Allmu1(:,:,i) = mu1;
    Allmu2(:,:,i) = mu2;
    if flag_calc_original_frangi == true
        S2_ori =Lambda1.^2 + Lambda2.^2;
        Ifiltered_ori =  exp(-Rb/beta).*(ones(size(I))-exp(-S2_ori/c));% 
        if(options.BlackWhite)
            Ifiltered_ori(Lambda2<0) = 0;
        else
            Ifiltered_ori(Lambda2>0) = 0;
        end
        ALLfiltered_ori(:,:,i) = Ifiltered_ori;
        if flag_use_original_frangi==true
            ALLfiltered(:,:,i) = Ifiltered_ori;
        end
        AllLambda1(:,:,i) = Lambda1;
        AllLambda2(:,:,i) = Lambda2;
    end
end

outDoH = Allmu1(:,:,1).*Allmu2(:,:,1);
for i = 2:length(sigmas)
    outDoH = max(outDoH, Allmu1(:,:,i).*Allmu2(:,:,i));
end

% Return for every pixel the value of the scale(sigma) with the maximum 
% output pixel value
if length(sigmas) > 1,
    [outIm,whatScale] = max(ALLfiltered,[],3);
    outIm = reshape(outIm,size(I));
    if(nargout>1)
        whatScale = reshape(whatScale,size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
    end
else
    outIm = reshape(ALLfiltered,size(I));
    if(nargout>1)
            whatScale = ones(size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles,size(I));
    end
end


