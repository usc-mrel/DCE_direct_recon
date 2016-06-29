% this demo is to calculate the backward modelling from image to TK
% parameter maps and from TK maps to image

% Yi Guo, 05/2016

%close all; clc; clear all;
load DCE50_0402;

[kx,ky,kz,nt,ncoil]=size(k);
X=sum(conj(repmat(sMaps,[1 1 1 nt 1])).*ifft2(k),5); % get fully-sampled DCE images
imgR=normalization(double(X));

ns=1; %get just one slice
imgR=abs(imgR(:,:,ns,:));
%% set time
delay=8; % delay frames for contrast injection
tpres=5/60; % temporal resolution, unit in minutes!

time=[zeros(1,delay),[1:(nt-delay)]*tpres];

alpha=pi*15/180; % flip angle, unit in radians
TR=0.006; % Repetition time, unit in seconds
Sb=imgR(:,:,:,1); % take pre-contrast first frame for baseline image
Sb=repmat(Sb,[1 1 1 nt]); % repmat in temporal dimension for substraction
%% calculate ktrans and test for/back-ward modelling
% Signal to concentration and to Ktrans
if ~exist('R1','var')  % use simulated uniform M0 and R1 if none exists
M0=5*ones(kx,ky,'single'); %use simulated M0, R1
R1=1*ones(kx,ky,'single');
end

CONC = sig2conc2(real(imgR),R1,M0,alpha,TR); % convert images to contrast concentration
[Kt,Vp]=conc2Ktrans(CONC,time); % convert concentration to TK maps
%% get back
% Ktrans back to concentration then to Signal, Sb is assumed known
[CR,A]=Ktrans2conc(Kt,Vp,time);
S = conc2sig(CR,R1,M0,Sb,alpha,TR);

% Use restored Signal S to do calculation again
CONC_r = sig2conc2(S,R1,M0,alpha,TR);
[Kt_r,Vp_r]=conc2Ktrans(CONC_r,time);

%% display
ns=1;
figure;imagesc(cat(1,cat(2,Kt(:,:,ns),Vp(:,:,ns)/2),cat(2,Kt_r(:,:,ns),Vp_r(:,:,ns)/2)),[0 0.1]);axis image; axis off; colorbar;

