function Pw = psdw(W,coi,f,dt)
%% PSDW Calculate power spectral density from wavelet coefficients
% 
%   Pw = PSDW(W,coi,f,dt,dj)
% 
%   INPUT:
%           W:   Wavelet coefficient matrix
%           coi: Cone of influence vector
%           f:   Wavelet frequency vector
%           dt:  Sampling period in units of s
% 
%   OUTPUT:
%           Pw: Power spectral density vector in units of X^2/Hz,
%               where X is unit of observed time series
% 
% Author: Michael von Papen
% 
% Date: 15.10.15

% Set sampling period delta t in units of s
if nargin<4
    dt = 1/2456;
end

% Convert wavelet coefficients to power
if ~isreal(W)
    W = 2*dt*abs(W).^2;
end

% Do not count coefficients within cone-of-influence
if nargin >= 3
    W = coi2nan(f,W,coi);
end

% Calculate PSD
Pw = nanmean(W,2);