function Cp=SAIF_p(tModel,A,T,sigma,alpha,beta,s,tau)
% function to generate parker's AIF 
% input tModel format: zeros in front to decide the length of delay
%                      tModel unit should be minutes

% Yi Guo, Yinghua Zhu, 10/25/2014

td=tModel>0;
if mean(diff(tModel))>0.5
tModel=tModel/60; % convert to minute unit
end

if nargin==1
A = [0.309 0.330]; %change A1 to 0.309 to decrease the peak, original [0.809 0.330]
T = [0.17046 0.365];
sigma = [0.0563 0.132];
alpha = 1.050;
beta = 0.1685;
s = 38.078;
tau = 0.483;
end

Cp = zeros(size(tModel));
for n = 1:2
    Cp = Cp + A(n)/sigma(n)/sqrt(2*pi)*exp(-(tModel-T(n)).^2/2/sigma(n)^2);
end
Cp = Cp + alpha*exp(-beta*tModel)./(1+exp(-s*(tModel-tau)));

Cp=Cp.*td;

Cp=Cp/3; % scale to mimic our data

end
