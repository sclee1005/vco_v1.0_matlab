function [Lambda1,Lambda2,Ix,Iy, mu1, mu2]=eig2imageModifiedBingxiongLin(Dxx,Dxy,Dyy)
% This function eig2image calculates the eigen values from the
% hessian matrix, sorted by abs value. And gives the direction
% along the ridge (eigenvector correspondes to the smallest eigenvalue) .
% 
% [Lambda1,Lambda2,Ix,Iy]=eig2image(Dxx,Dxy,Dyy)
% small lamda correspondes to vectors along vessel, large lamda across
% vessel
% abs(Lambda1)<abs(Lambda2), mu1<mu2
%
% | Dxx  Dxy |
% |          |
% | Dxy  Dyy |


% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v1x = 2*Dxy; v1y = Dyy - Dxx + tmp;

% Normalize
mag = sqrt(v1x.^2 + v1y.^2); i = (mag ~= 0);
v1x(i) = v1x(i)./mag(i);
v1y(i) = v1y(i)./mag(i);

% The eigenvectors are orthogonal (only for symmetric matrix)
v2x = -v1y; 
v2y = v1x;

% Compute the eigenvalues
mu2 = 0.5*(Dxx + Dyy + tmp); % mu1<mu2 
mu1 = 0.5*(Dxx + Dyy - tmp); % negative: bright tubular:   positive: dark tubular

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);

Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);

Ix=v2x; Ix(check)=v2x(check); %this is related to Hessian2D.m  [Y,X] = ndgrid
Iy=v2y; Iy(check)=v2y(check);






