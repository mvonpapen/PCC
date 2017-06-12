%% Function to compute PLI from wavelet cross spectrum

function PLI = phase_lag_index ( Im, scale, nsig, dt, varargin )

% input imaginary part of coherency matrix for wPLI [Vinck et al., 2011]
args = struct('wPLI', false);
        
args = parseArgs(varargin,args);
wPLI = args.wPLI;

if nargin<4
    dt = 1/2456;
end
if nargin<3
    nsig = 6;
end

% Average sign of phases over time
if ~wPLI
    Ph_avg = temp_avg_wave( sign(Im), dt, nsig*scale );
else
    % biased est
    Ph_avg_num   = temp_avg_wave( Im, dt, nsig*scale );
    Ph_avg_denom = temp_avg_wave( abs(Im), dt, nsig*scale );
    Ph_avg       = Ph_avg_num./Ph_avg_denom;
%     % debiased est (taken from fieldtrip
%     Im      = imag(Wxy);
%     outsum  = temp_avg_wave( Im, dt, nsig*scale );
%     outsumW = temp_avg_wave( abs(Im), dt, nsig*scale );
%     outssq  = temp_avg_wave( Im.^2, dt, nsig*scale, 'smooth','gauss_square' );
%     Ph_avg  = (outsum.^2 - outssq)./(outsumW.^2 - outssq);
end
  
PLI = abs(Ph_avg);