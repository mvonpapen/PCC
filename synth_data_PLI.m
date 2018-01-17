%% Parameter
nsig  = 6;
w0    = 12;
sig   = 0.705; %p01=0.705, p05=0.555
N     = 10;
% NOTE: If wPLI=true use Im(Cohy), else Im(Wxy)
wPLI  = false;

% Parameters for synth data
alpha = 10;  %coh
beta  = 20;  %vc
gamma = 50;  %coh
f     = 1:70;
dt    = 1/2500;
nt    = 20/dt;
t     = (0:nt-1)*dt;
nl    = 3;
dp    = 30/180*pi;
A     = 1000; % amplitude correction

% Filename to save data to
fnam  = ['PSD_synth_PLI_ns' mat2str(nsig) '_w' mat2str(w0) '_nl' mat2str(nl) '.mat'];

%% Preallocate
Ptot = zeros(length(f),N);
Pinc = zeros(length(f),N);
Pcoh = zeros(length(f),N);

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
    Wxy              = squeeze(W1(:,:,1).*conj(W1(:,:,2)));
    PLI              = phase_lag_index(imag(Wxy), scale, nsig, dt, 'wPLI', wPLI);
    [a, b]           = pcc ( f, W1, repmat(PLI,1,1,2,2), coi, sig );
    Pcoh(:,j)        = a(:,1,2);
    Pinc(:,j)        = b(:,1,2);
    
    fprintf('\b\b\b\b\b\b\b\b\b')
end

% % Save if desired
% fprintf('\n Saving results in \n %s\n', [pwd '\' fnam])
% save(fnam, 'A', 'N', 'sig', ...
%     'Pcoh', 'Pinc', 'Ptot', 'alpha', 'beta', 'gamma', 'nt', ...
%     'coi', 'f', 'w0', 'nsig', 'nl', 'dp')
