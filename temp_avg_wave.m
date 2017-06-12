function avg_wave = temp_avg_wave(wave,dt,speriod,varargin)
%% TEMP_AVG_WAVE Function used to average wavelet matrix over time, input parameters are:
% 
%   INPUT:
%           wave:   (matrix from wavelet transform)
%           dt:     (time step)
%           speriod:(vector of periods used for smoothing, e.g. speriod=1./f. 
%                   Must have the same number of elements as the frequency 
%                   vector used for wavelet transformation. speriod is in
%                   unit of seconds!)
%
% 	OPTIONAL INPUT:
%           'kernel':       ('gauss' or 'box', default 'gauss', type of 
%                           smoothing window used for smoothing in time)
% 
%   OUTPUT:
%           avg_wave: Averaged wavelet coefficient matrix
% 
% Author: Michael von Papen
% 
% Date: 14.10.15


%% parameters
speriod=speriod(:);
[m, n] = size(wave);
avg_wave = zeros(m,n);

args = struct('kernel','gauss');
args = parseArgs(varargin,args);


%% creating frequency vector for fft
f = [0:fix(n/2), -(fix((n-1)/2):-1:1)]'/(n*dt);

%% Take out zero frequency so that ifft does not produce imaginary values
tag = 0;
if speriod(1) == Inf
    wave0   = wave(1,:);
    wave    = wave(2:end,:);
    speriod = speriod(2:end);
    tag     = 1;
end


%% actual smoothing for each type by convolution using fft
switch upper(args.kernel)
    case 'GAUSS'
        F = exp( - 2 * pi^2 * f.^2 * speriod'.^2  ); % = fft of normalized gaussian with std = speriod
        avg_wave = ifft(F.*fft(wave'))';
    case 'BOX'
        F = sinc( f * speriod' ); % = fft of normalized boxcar with width = speriod
        avg_wave = ifft(F.*fft(real(wave')))'+1i*ifft(F.*fft(imag(wave')))';
end


%% Insert zero frequency value
if tag == 1
    avg_wave = [wave0; avg_wave];
end