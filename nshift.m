function k=nshift(k,dim)
% fftshift the specific dimension

% Yi Guo 08/14

for i=1:length(dim)
    k=fftshift(k,dim(i));
end
end