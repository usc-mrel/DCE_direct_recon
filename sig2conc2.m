function C = sig2conc2(img,R10,M0,alpha,TR,imgB)
% functon to calculate contrast concentration from image intensity
% equation: 
% R1(t)=-1/TR*ln(1-((S(t)-S(0))/S0*sin(alpha))+(1-m)/(1-m*cos(alpha)))
%over 1-cos(alpha)*((S(t)-S(0))/S0*sin(alpha))+(1-m)/(1-m*sin(alpha)))
% where m=exp(-R10*TR)
% then C(t)=(R1(t)-R1(0))/r1

% Yi Guo, 06/12/2014
% No mask comparing to Marc's version, otherwise allmost the same
% some simplification, pay attention to R1, and alpha unit!!!

if nargin==5
    imgB=img(:,:,:,1);
end

Rcs=4.39; % contrast relaxivity
nt=size(img,4); % temporal dimension size

m=exp(-repmat(R10,[1 1 1 nt])*TR);
par2=(1-m)./(1-m*cos(alpha));

par1=(img-repmat(imgB,[1 1 1 nt]))./repmat(M0+eps,[1 1 1 nt])/sin(alpha);

B=(par1+par2);

Rt=-1/TR*real(log((1-B)./(1-cos(alpha)*B+eps)));

R1B=Rt(:,:,:,1); % baseline R1 is equal to R10 if img(1)==imgB
%R1B=R10;
C=Rt-repmat(R1B,[1 1 1 nt]);

C=C/Rcs;

C(C<0) = 0;  % get rid of outliers
C(C>180)=0;
end

    
    



