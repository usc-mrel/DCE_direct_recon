function [Kt,Vp,iter,fres]=Kt_Vp_SEN(Kt,Vp,sMaps,U1,kU,opt)
% This is the function to alternatively reconstruct Ktrans and Vp using a
% small number of l-BFGS iterations
% Input: Kt: Initial Ktrans guess
%        Vp: Initial Vp guess
%        sMaps: sensitivity maps
%        U1: under-sampling mask
%        kU: under-sampled k-space
%        opt: other parameters
%Output: Kt,Vp, Kt_inter,Vp_inter: reconstructed and intermediate TK maps
%        fres: objective function value across iterations
%        gradres: gradient value across iterations

% Yi Guo, 06/2014

options.MaxIter=opt.Initer;
options.display = 'off';
options.Method ='lbfgs'; % use l-BFGS method in minFunc
%options.optTol=1e-3;
%options.progTol=5e-4;
options.useMex=0;
options.inter=0;
options.numDiff=0;

exitflag=0;
iter=1;
fres=[];

while(exitflag==0)

opt.Vp=reshape(Vp,[opt.size(1) opt.size(2)]);
opt.lambda1=opt.lambdaA(1);
opt.lambda2=opt.lambdaA(2);
[Kt,f,exitflag,output]=minFunc(@Ktrans2sig_sen_WT,Kt(:),options,sMaps,U1,kU,opt);

fres=[fres;output.trace.fval];

opt.Kt=reshape(Kt,[opt.size(1) opt.size(2)]);
opt.lambda1=opt.lambdaA(3);
opt.lambda2=opt.lambdaA(4);

[Vp,f,exitflag,output]=minFunc(@Vp2sig_sen_WT,Vp(:),options,sMaps,U1,kU,opt);

fres=[fres;output.trace.fval];

if opt.plot
    imagesc(real(cat(2,opt.Kt,opt.Vp/2)),[0 0.1]);axis image; axis off; 
    title(['Ktrans, Vp/2, Out iter=',num2str(iter)]); colorbar; 
    drawnow;
end

iter=iter+1;
end

Kt=reshape(Kt,[opt.size(1) opt.size(2)]);
Vp=reshape(Vp,[opt.size(1) opt.size(2)]);


end

