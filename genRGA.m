function [petable, kspace] = genRGA(FOVy, FOVz, ky, kz, totlength, flag_option, precontrast, splitN, randprob, pertangle, hdRESy, period, randseed, pertseed, initangle)

% [petable, kspace] = genRGA(FOVy, FOVz, ky, kz, totlength, flag_option, precontrast, splitN, randprob, pertangle, hdRESy, period, randseed, pertseed, initangle)
% Generate Cartesian-radial Randomized Golden Angle phase encode table
% 	
% Output the table as a file.
% 	
% Inputs:
% 	FOVy:			(default 240.0) FOV along y in mm
% 	FOVz:			(default 240.0) FOV along z in mm
% 	ky:				(default 256) dimension of matrix in y
% 	kz:				(default 128) dimension of matrix in z
% 	totlength:		(default 50000) total number of output phase encodes
% 	flag_option:	(default '0000') binary code for all mode combinations
% 		-(0001) spokes reaches the boundary of an ellipse bounded by the matrix
%		-(0010) spoke out-in to reduce eddy current
%		-(0100) output a separate file of sampling table
%		-(1000) flip every second phase encode
%   precontrast:    (default 0.0) [0,1], portion of full-frame data to acquire as reference
%	splitN:         (default 1) the number of fold for output kspace splitting 
% 	randprob:		(default 1.0) threshold of the probability to sample a phase encode
%	pertangle:      (default 0.0) perburb angle in radius/pi, will perturb half angle each side
% 	hdRESy:			(default 0.0) high-density resolution along y in mm
% 	period:			(default 0) the number of TR's within which a PE is not revisted
% 	randseed:		(default 1) seed for generating a pseudorandom series
%	pertseed:       (default 1) seed for generating a pseudorandom angle
% 	initangle:		(default 0.0) initial angle of the radial in radians
%
% Outputs:
%   petable:  sampling petable (ky, kz)
%   k_space:  sampling map ([ky,kz] without split, [ky,kz,splitN] with split w/o precon, [ky,kz,splitN+1] with split and precon)
%   and a separate file including output_z and output_y
%
% Example:
%   [petable, kspace1] = genRGA(240.0, 240.0, 256, 186, 50000, bin2dec('1111'), 0.25, 5, 0.3, 0.1, 10.0, 2500, 1, 1, 0.0);
% or
%   [petable, kspace2] = genRGA;
%
% Use
%   for n=1:6, subplot(2,3,n), imshow(kspace1(:,:,n)); end
% or
%   imshow(kspace2);
% to see the sampling maps.
%
% 	--------------------
% 	Version 1.1, 07/21/2014
%
% 	Authors: Yinghua Zhu (yinghuaz@usc.edu)
%            Yi Guo
%            Marc Lebel
% 			 Krishna Nayak

% parameter checking
% A. set default input values
if (nargin < 1),    FOVy = 240.0;       end
if (nargin < 2),    FOVz = 240.0;       end
if (nargin < 3),    ky = 256;           end
if (nargin < 4),    kz = 186;           end
if (nargin < 5),    totlength = round(ky*kz/5);     end
if (nargin < 6),    flag_option = bin2dec('0001');  end
if (nargin < 7),    precontrast = 0.0;  end
if (nargin < 8),    splitN = 1;         end
if (nargin < 9),    randprob = 0.3;     end
if (nargin < 10),   pertangle = 0.0;    end
if (nargin < 11)
    R = (ky*kz) / round((totlength-ky*kz*precontrast)/splitN);
    hdRESy = FOVy/ky*sqrt(R*10);
end
if (nargin < 12),   period = round(ky*kz/50*0.8);   end
if (nargin < 13),   randseed = 1;       end
if (nargin < 14),   pertseed = 1;       end
if (nargin < 15),   initangle = 0.0;    end

% B. safety check - to be more checks
ky = floor(ky/2)*2;
kz = floor(kz/2)*2;
period = min(period, ky*kz/2);
if precontrast>1,   precontrast = 1;    end
if randseed==0,   randseed = 1;         end
if pertseed==0,   pertseed = 1;         end

% C. funtion parameter define
flag_ellipse    = bitand(flag_option, 1) > 0;
flag_out2in     = bitand(flag_option, 2) > 0;
flag_write      = bitand(flag_option, 4) > 0;
flag_pair       = bitand(flag_option, 8) > 0;
flag_precon     = (precontrast ~= 0);
flag_split      = (splitN > 1);
flag_probability = (randprob < 1);
flag_perturb    = (pertangle > 0);
flag_center     = (hdRESy > 0);
flag_efficient  = (period > 0);
% done input parameter checking



% initialization
petable = zeros(totlength, 2);
kspace = zeros(ky, kz);
counter = 1;            % output counter
RESy = FOVy/ky;
RESz = FOVz/kz;
centerrange = RESy/hdRESy;
yzFOVratio = FOVy/FOVz;
yzRESratio = RESy/RESz;
randrange = 10000;

% golden angle
ga = pi*(3-sqrt(5));    % golden angle
rlength = (ky+kz)/2;    % long enough radial
nradial = 0;
if flag_pair
    ga = ga/2;
end

% initialization, fixed length
sample_index = zeros(1, rlength);
peflag_rectangle = ones(1, rlength);
peflag_ellipse = ones(1, rlength);
peflag_probability = ones(1, rlength);
peflag_center = zeros(1, rlength);

if flag_ellipse
    c = sqrt(abs( (kz/2)^2 - (ky/2)^2 ));       % for ellipse
end

if flag_efficient
    sample_record = ones(ky, kz)*(totlength + period);  % initialization
end

if flag_probability
    newrandseed = myrand(randseed);
end

if flag_perturb
    newpertseed = myrand(pertseed);
end

if flag_precon
    yPrecon = floor(floor(sqrt(precontrast)*ky)/2)*2;
    zPrecon = floor(floor(sqrt(precontrast)*kz)/2)*2;
    offset = yPrecon*zPrecon;
    for m = 1:zPrecon
        for n = 1:yPrecon
            if totlength>0
                petable((m-1)*yPrecon+n, 1) = n + (ky-yPrecon)/2;
                petable((m-1)*yPrecon+n, 2) = m + (kz-zPrecon)/2;
                kspace(petable((m-1)*yPrecon+n,1), petable((m-1)*yPrecon+n,2)) = kspace(petable((m-1)*yPrecon+n,1), petable((m-1)*yPrecon+n,2)) + 1;
                totlength = totlength-1;
            end
        end
    end
else
    offset = 0;
end

if flag_split
    lengthN = floor(totlength/splitN);
    indN = 1;       % index for # of split
    cntN = 0;       % counter inside each split
    if flag_precon
        kspaceN = zeros(ky, kz, splitN+1);
        kspaceN(:,:,indN) = kspace;
        indN = indN + 1;
        kspace = zeros(ky, kz);     % reset
    else
        kspaceN = zeros(ky, kz, splitN);
    end
end 

% generating spoke
while totlength>0
    radial_ga = rlength * exp(1i*(ga*nradial+initangle));
    angle_ga = angle(radial_ga);
    
    % angle adjustment for different FOV case
    if yzFOVratio ~= 1
        angle2_ga = (atan(tan(angle_ga) * yzFOVratio));     % adjusted angle in diffFOV case
        if cos(angle_ga)<0
            angle2_ga = angle2_ga + pi;
        end
        radial2_ga = rlength * exp(1i* angle2_ga);
    else
        radial2_ga = radial_ga;
        angle2_ga = angle_ga;
    end
    
    % core: Bresenham's line search algorithm (round to N+0.5)
    [z, y] = bresenham(-sign(real(radial2_ga))*.01+.5, -sign(imag(radial2_ga))*.01+.5, real(radial2_ga)+.5, imag(radial2_ga)+.5);
    % z, y are one entry longer on both ends
    z = z(2:end-1);
    y = y(2:end-1);
    radial_length = length(z);
    
    % perturb a radial randomly
    if flag_perturb
        for n_perturb = 1:length(z)
            perb = z(n_perturb) + 1i*y(n_perturb);
            perb = perb * exp( (mod(newpertseed,randrange)/randrange - 1) * pertangle*1i*pi * cos(angle2_ga) / cos(angle_ga) );
%             perb = perb * exp( (mod(newpertseed,randrange)/randrange - 1) * pertangle*1i*pi * cos(angle2_ga) / cos(angle_ga) * yzFOVratio); 
            y(n_perturb) = round(imag(perb));
            z(n_perturb) = round(real(perb));
            newpertseed = myrand(newpertseed);
        end
    end
    
    z = z - .5;
    y = y - .5;
    
    if flag_ellipse         % ellipse mode
        if ky > kz
            dist = sqrt(z.^2 + (y-c).^2) + sqrt(z.^2 + (y+c).^2);
            peflag_ellipse(dist > ky) = 0;
            peflag_ellipse(dist <= ky) = 1;
        else
            dist = sqrt((z-c).^2 + y.^2) + sqrt((z+c).^2 + y.^2);
            peflag_ellipse(dist > kz) = 0;
            peflag_ellipse(dist <= kz) = 1;
        end
        temp = find(peflag_ellipse(1:radial_length) == 1);
        radial_length = temp(end);
    else                    % regular rectangle support
        peflag_rectangle = zeros(1, rlength);
        peflag_rectangle(abs(z) <= kz/2-0.5) = 1;
        peflag_rectangle(abs(y) > ky/2-0.5) = 0;
        temp = find(peflag_rectangle(1:radial_length) == 1);
        radial_length = temp(end);
    end
    
    % probability mode
    if flag_probability
        for m = 1:radial_length
            if (mod(newrandseed,randrange)/randrange) >= randprob
                peflag_probability(m) = 0;
            else
                peflag_probability(m) = 1;
            end
            newrandseed = myrand(newrandseed);
        end
    end
    
    % center mode
    if flag_center
        % consider yzRESratio, similar to processing after bresenham
        y1=y; z1=z;
        if 1    % center uses ellipse anyway
            c1 = sqrt(abs( (kz*centerrange/yzRESratio/2)^2 - (ky*centerrange/2)^2 ));
            if ky > kz/yzRESratio
                dist = sqrt(z1.^2 + (y1-c1).^2) + sqrt(z1.^2 + (y1+c1).^2);
                peflag_center(dist <= ky*centerrange) = 1;
                peflag_center(dist > ky*centerrange) = 0;
            else
                dist = sqrt((z1-c1).^2 + y1.^2) + sqrt((z1+c1).^2 + y1.^2);
                peflag_center(dist <= kz*centerrange/yzRESratio) = 1;
                peflag_center(dist > kz*centerrange/yzRESratio) = 0;
            end
        else
            peflag_center = ones(1, rlength);
            peflag_center(abs(z1) > kz/2*centerrange/yzRESratio-0.5) = 0;
            peflag_center(abs(y1) > ky/2*centerrange-0.5) = 0;
        end
    end
    
    % consider for efficient mode
    n = 1;
    for m = 1:radial_length
        if ( (peflag_rectangle(m) && peflag_ellipse(m)) && (peflag_probability(m) || peflag_center(m)) )
            if flag_efficient
                if (sample_record( ky/2+0.5-y(m), kz/2+0.5+z(m) ) - (totlength - n)  >= period)
                    % if bigger than the gap, sample it !
					sample_index(n) = m;
					% update the sample_record
					sample_record( ky/2+0.5-y(m), kz/2+0.5+z(m) ) = totlength - n;
					n = n+1;
                end
            else
                sample_index(n) = m;
                n = n+1;
            end
        end
    end
    
    % update the length
    radial_length = n-1;
    
    if radial_length > 0
        temp = totlength;
        totlength = totlength - radial_length;
        if (totlength<0)
            radial_length = temp;
        end
        
        if flag_out2in
            outputOrder = radial_length:-1:1;   % outside in
        else
            outputOrder = 1:radial_length;      % inside out
        end

        for m = outputOrder
            petable(counter+offset, 1) = ky/2+0.5-y(sample_index(m));
            petable(counter+offset, 2) = kz/2+0.5+z(sample_index(m));
            
            if flag_pair
                if mod(counter,2) == 1
                    petable(counter+offset, 1) = ky+1-petable(counter+offset, 1);
                    petable(counter+offset, 2) = kz+1-petable(counter+offset, 2);
                end
            end
            
            kspace(petable(counter+offset,1), petable(counter+offset,2)) = kspace(petable(counter+offset,1), petable(counter+offset,2)) + 1;
            counter = counter +1;
            
            % consider for splitN
            if flag_split
                cntN = cntN + 1;
                if cntN == lengthN
                    kspaceN(:,:,indN) = kspace;
                    indN = indN + 1;
                    cntN = 0;                 % reset
                    kspace = zeros(ky, kz);     % reset
                end
            end
        end
    end
    
    nradial = nradial+1;
end

% if split, use splited kspace as output
if flag_split
    kspace = kspaceN;
end

if flag_write
    dlmwrite('RGApetable.txt',petable, ' ');
end

end
% genRGA ends here


function randno = myrand(seed)
% My pseudorandom number generator, which was used in Matlab.
% Return the next pseudorandom given a seed.
% 
% Reference:
% http://www.mathworks.com/moler/random.pdf
% S. K. Park and K. W. Miller, Random number generators: Good ones are hard
% to find, Communications of the ACM, 31 (1988), pp. 1192-1201.
% ZYH, 02/25/2014
next = seed * 7^5;
randno = mod(next, 2^31-1);
end

function [x, y] = bresenham(x1, y1, x2, y2)

% function [x, y] = bresenham(x1, y1, x2, y2)
% Bresenham line search algorithm
%
% Inputs:
%   (x1,y1): Start position
%   (x2,y2): End position
%
% Outputs:
%   (x, y): the line coordinates from (x1,y1) to (x2,y2)
%
% Modified from online code
% Yinghua Zhu, 08/16/2013

% If not steep, change to non-steep
dx = abs(x2-x1);
dy = abs(y2-y1);
steep = dy > dx;
if steep
    dy = dx;
end

% The main algorithm goes here
% Calculate the length, minimum is 1
if steep
    qlength = abs(round(y1)-round(y2))+1;
else
    qlength = abs(round(x1)-round(x2))+1;
end

if dy == 0
    q = zeros(qlength, 1);
else
    q = (0: dy/(qlength-1): dy)';
end

% calculate coordinates for each point
if steep
    if y1 <= y2
        y = (round(y1):round(y2))';
    else
        y = (round(y1):-1:round(y2))';
    end
    if x1 <= x2
        x = x1+q;
    else
        x = x1-q;
    end
else
    if x1 <= x2
        x = (round(x1):round(x2))';
    else
        x = (round(x1):-1:round(x2))';
    end
    if y1 <= y2
        y = y1+q;
    else
        y = y1-q;
    end
end

% round to grid points for outputs
x = round(x);
y = round(y);
end