%% Parameter
nsig  = 6;
w0    = 12;
sig   = sig_coh_thresh(w0, nsig);
phase_thresh = 15.5;
N     = 10;

% Parameters for synth data
alpha = 10;  %coh
beta  = 20;  %vc
gamma = 50;  %coh
dt    = 1/2500;
T     = 20; %length of time series in s
f     = 1:100;
% f     = logspace(log10(1/T),log10(1/2/dt),400); # for reconstruction
nt    = T/dt;
t     = (0:nt-1)*dt;
nl    = 3;
dp    = 30/180*pi;
A     = 1000; % amplitude correction

% Filename to save data to
fnam  = ['PSD_synth_PCC_ns' mat2str(nsig) '_w' mat2str(w0) '_nl' mat2str(nl) '.mat'];


%% Preallocate
Ptot = zeros(length(f),N);
Pinc = zeros(length(f),N);
Pcoh = zeros(length(f),N);
Pvc  = zeros(length(f),N);


fprintf ('\r Processing ')
for j=1:N
    fprintf('%4d/%4d', j, N)
    
    %% Synthetic time series
    [x, y] = synth_data( t, nl, A, dp );


    %% Spectral analysis
    scale            = (w0+sqrt(2+w0^2))/4/pi ./ f;
    [X, W1, coi, P]  = preprocdata([x y], 'freq', f, 'w0', w0, 'dt', dt);
    sigma2           = var(X(:,1)); % variance of original time series
    Ptot(:,j)        = P(:,1);
    [C, Wxy, W]      = wave_coherence(W1, scale, nsig, 1, dt);
    [a, b, c]        = pcc ( f, W, C, coi, sig, Wxy, phase_thresh );
    Pcoh(:,j)        = a(:,1,2);
    Pinc(:,j)        = b(:,1,2);
    Pvc(:,j)         = c(:,1,2);
    
    fprintf('\b\b\b\b\b\b\b\b\b')

end

% % Save if desired
% fprintf('\n Saving results in \n %s\n', [pwd '\' fnam])
% save(fnam, 'A', 'N', 'phase_thresh', 'w0', 'nsig', 'nl', 'dp', 'sig', ...
%    'Pcoh', 'Pinc', 'Pvc', 'Ptot', 'alpha', 'beta', 'gamma', ...
%    'coi', 'f')