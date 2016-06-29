%% load data
close all; 
clc; clear all;

addpath(genpath('minFunc_2012'));

load DCE50_0402
%load T1_0402;

ns=1; % choose one slice of k-space
k=k(:,:,ns,:,:);
opt.size=size(k);  
[kx,ky,kz,nt,ncoil]=size(k);

sMaps=sMaps(:,:,ns,:,:);
sMaps=reshape(sMaps,[kx ky 1 1 ncoil]);

if ~exist('R1','var')  % use simulated uniform M0 and R1 if none exists
M0=5*ones(kx,ky,'single'); %use simulated M0, R1
R1=1*ones(kx,ky,'single');
end

imgF=sum(conj(repmat(sMaps,[1 1 1 nt 1])).*ifft2(k),5); % get fully-sampled

%% set parameters
opt.wname='db4'; % wavelet parameters
opt.worder={[1 2],[1 2],[1 2]};

opt.R1=R1;
opt.M0=M0;
opt.Sb=repmat(imgF(:,:,:,1),[1 1 1 nt]);  %baseline image
opt.alpha=pi*15/180; %flip angle
opt.TR=0.006;  %TR

delay=8; % delay frames for contrast injection
tpres=5/60; % temporal resolution, unit in seconds!
opt.time=[zeros(1,delay),[1:(nt-delay)]*tpres];
opt.plot=1;  % to plot intermediate images during recon

opt.lambdaA=[0.000 0.000 0.000 0.000]; % Kt:TV, Wavelet, Vp: TV, wavelet
opt.Initer=10;  %inter interations
opt.Outiter=10; % outer interations
%% calculate fully-sampled Ktrans and Vp
CONCF = sig2conc2(real(imgF),R1,M0,opt.alpha,opt.TR);
opt.AIF=SAIF_p(opt.time); % get population-averaged AIF
[Kt,Vp]=conc2Ktrans(CONCF,opt.time,opt.AIF);
%% undersamping by RGR
R = 30; %under-sampling rate
[~, U11] = genRGA(220, 220, kx,ky, round(kx*ky/R*nt), bin2dec('0001'), 0, nt);

U1=reshape(nshift(U11>0,[1 2]),[kx,ky,1,nt]);
U1=repmat(U1,[1 1 1 1 ncoil]);
U1(:,:,:,1,:)=1;  % first frame is fully-sampled
kU=k.*U1;
imgU=sum(conj(repmat(sMaps,[1 1 1 nt 1])).*ifft2(kU),5);
CONC1 = sig2conc2(real(imgU),R1,M0,opt.alpha,opt.TR);
[Kt_U,Vp_U]=conc2Ktrans(CONC1,opt.time,opt.AIF);
%% Direct reconstruction
% use zeros as initial guess
Kt_1=zeros(opt.size(1),opt.size(2)); 
Vp_1=Kt_1;

tic,
[Kt_r1,Vp_r1,iter,fres]=Kt_Vp_SEN(Kt_1,Vp_1,sMaps,U1,kU,opt);
toc,

iter,
Kt_r1=real(Kt_r1);  % get real value after recon
Vp_r1=real(Vp_r1);

%% display results
figure;
subplot(2,3,1); imagesc(Kt,[0 0.1]); title('Fully-sampled Ktrans'); colorbar;axis image; axis off;
subplot(2,3,2); imagesc(Kt_U,[0 0.1]); title('Zero-padded Ktrans'); colorbar;axis image; axis off;
subplot(2,3,3); imagesc(Kt_r1,[0 0.1]); title('Reconstructed Ktrans'); colorbar;axis image; axis off;
subplot(2,3,4); imagesc(Vp,[0 0.2]); title('Fully-sampled Vp'); colorbar;axis image; axis off;
subplot(2,3,5); imagesc(Vp_U,[0 0.2]); title('Zero-padded Vp'); colorbar;axis image; axis off;
subplot(2,3,6); imagesc(Vp_r1,[0 0.2]); title('Reconstructed Vp');  colorbar;colorbar;axis image; axis off;

figure;
plot(fres);title('Objective function changes over iteration');
xlabel('Iteration'); ylabel('Objective function value');

