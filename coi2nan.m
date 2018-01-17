function waveOut=coi2nan(f,waveIn,coi)

%%COI2NAN Set data within the cone-of-influence (coi) to NaN
% waveOut=coi2nan(f,waveIn,coi)
%
%   IN: 
%       f:      frequency vector
%       waveIn: Wavelet coefficient Matrix
%       coi:    cone of influence: vector with minimal 
%               periods that are not affected by border
%   OUT:
%       waveOut: Wavelet coefficient Matrix with NaNs
%
%   Author: Michael von Papen, Date: 13.10.15

[a, b, c] = size( waveIn );

if nargin < 3; coi = gencoi(b); end

% Check size of coi
[b2, c2] = size( coi );
if b ~= b2
    coi = coi';
    [b2 c2] = size( coi );
    if b ~= b2
        error('Number of data points don''t match!')
    end
    if c2~=1
        error('COI must be the same for all channels.')
    end
end

% Generate coi matrix
C = repmat(coi(:)',[a, 1, c]);

% Generate matrix of periods
warning('off') % switch off warnings to allow divide by zero for f=0
P = repmat(1./f(:),[1, b, c]);
warning('on')

% Use coi matrix and periods as mask for wavelet coefficients
waveOut      = waveIn;
waveOut(C<P) = NaN;