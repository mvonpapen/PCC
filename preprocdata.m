function [ data, W, coi, Pw ] = preprocdata( data, varargin )
%% PREPROCDATA Preprocess data for PCC:
%  1. Compute wavelet transform
%  2. Compute PSD of both methods
%
%   INPUT:
%           data:    signal of k channels (nt x nch elements)
%
%   OPTIONAL INPUT:
%           dt:     Time increment
%           freq:   Frequency vector in Hz (nf elements) that shall be 
%                   evaluated
%
%  	OUTPUT:
%           W:      Wavelet transform of data (nf x nt x nch elements)
%           coi:    Wavelet cone of influence
%           Pw:     Global wavelet spectrum (nf x nch elements)
% 
% Author: Michael von Papen
% 
% Date: 14.10.15


args = struct('freq',logspace(0,2.7,30),...
            'dt',1/2456,...
            'w0', 12);
args = parseArgs(varargin,args);
arts = args.art;
w0 = args.w0;
dec = args.dec;


%% Parameters
[n, nch] = size(data);
if nch>n
    data = data';
    [n, nch] = size(data);
end
f = args.freq;
m = length(f);
dt = args.dt;


%% Preset variables
if nargout > 1
    W = NaN(m,n,nch);
    Pw = NaN(m,nch);
end


% Wavelet Transform
if nargout > 1
    W = NaN(m,n,nch);
    Pw = NaN(m,nch);
    coi = NaN(n,nch);
    for i = 1:nch
        [W(:,:,i),~,~,coi(:,i)] = waveTF_TC98(data(:,i),dt,f,1,'MORLET',w0);
        Pw(:,i) = psdw(squeeze(W(:,:,i)),coi(:,i),f);
    end
end