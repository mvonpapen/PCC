function PLI = phase_lag_index ( Im, scale, nsig, dt, varargin )

%%PHASE_LAG_INDEX Computes time-frequency resolved phase lag index
% PLI = phase_lag_index ( Im, scale, nsig, dt, varargin )
%
%   IN: 
%       IM:     imaginary part of wavelet cross spectrum
%       scale:  wavelet scale according to frequency f
%       nsig:   n_sigma, number of Gaussian window std used for averaging
%       dt:     Sampling period in units of s
%       varargin: 'wPLI', true or false. 
%                 If true, weighted PLI is calculated, otherwise the
%                 normal PLI is calculated
%   OUT:
%       PLI:    time-frequency resolved (weighted) Phase Lag Index [Nf x Nt]
%
%   Author: Michael von Papen, Date: 12.06.17
args = struct('wPLI', false);


% Parse argument (wPLI or PLI)
args = parseArgs(varargin,args);
wPLI = args.wPLI;


% Set sampling period
if nargin<4
    dt = 1/2456;
end


% Set n_sigma
if nargin<3
    nsig = 6;
end


%% Calculate (weighted) PLI
% Average sign of phases over time
if ~wPLI
    Ph_avg = temp_avg_wave( sign(Im), dt, nsig*scale );
else
    % biased est
    Ph_avg_num   = temp_avg_wave( Im, dt, nsig*scale );
    Ph_avg_denom = temp_avg_wave( abs(Im), dt, nsig*scale );
    Ph_avg       = Ph_avg_num./Ph_avg_denom;
%     % debiased est (taken from fieldtrip)
%     Im      = imag(Wxy);
%     outsum  = temp_avg_wave( Im, dt, nsig*scale );
%     outsumW = temp_avg_wave( abs(Im), dt, nsig*scale );
%     outssq  = temp_avg_wave( Im.^2, dt, nsig*scale, 'smooth','gauss_square' );
%     Ph_avg  = (outsum.^2 - outssq)./(outsumW.^2 - outssq);
end
PLI = abs(Ph_avg);