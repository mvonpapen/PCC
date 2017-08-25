%% Plot the outcome of several trials trying to replicate PSD of synthetic time series


%%  Choose a method (PCC, ImCoh, wPLI)
% % PCC
% meth = 'PCC';
% synth_data_PCC %generates synthetic data and applied PCC

% PLI
meth = 'wPLI';
synth_data_PLI %generates synthetic data and applied PCC based on wPLI



%% Figure
fig1 = figure('Papersize', [16 10], 'PaperPosition', [0.75 0.5 14.5 9], ...
        'PaperPositionmode', 'manual', 'Visible', 'on'); 
    
%% Parameters
ds    = 1;
ds2   = 150; %downsampling for plot
dt    = 1/2500;
nt    = 20/dt;
t     = (0:nt-1)*dt;
t     = t(1:ds:end);
limy  = [1e-8 1e-6];


%% Compute mean and error
% PCC
if strcmp(meth,'PCC')
    M = [mean(Ptot,2) mean(Pinc,2) mean(Pcoh,2) mean(Pvc,2)]';
    lineprops.col={'k'; 'r'; 'b'; 'g'};
    E = [ std(Ptot,0,2) std(Pinc,0,2) std(Pcoh,0,2) std(Pvc, 0,2) ]';
end
% IC or PLI
if any(strcmp(meth,{'wPLI', 'PLI', 'IC'}))
    M = [mean(Ptot,2) mean(Pinc,2) mean(Pcoh,2)]';
    lineprops.col={'k'; 'r'; 'b'};
    E = [ std(Ptot,0,2) std(Pinc,0,2) std(Pcoh,0,2) ]';
end

%% Plot results
mseb(f, M, E, lineprops);
xlim([1 60])
hold all
xa    = [zeros(nt/5,1); ...
         sin(2*pi*alpha*(1:3*nt/5)'*dt); zeros(nt/5,1)];
xb1   = [sin(2*pi*beta*(1+(1:nt/2)'/nt*0.25).*(1:nt/2)'*dt); ...
         sin(2*pi*beta*(1+(nt/2+1:nt)'/nt*0.25).*(nt/2+1:nt)'*dt-dp)*0];
xb2   = [sin(2*pi*beta*(1+(1:nt/2)'/nt*0.25).*(1:nt/2)'*dt)*0; ...
         sin(2*pi*beta*(1+(nt/2+1:nt)'/nt*0.25).*(nt/2+1:nt)'*dt-dp)];
xg    = [sin(2*pi*gamma*t(1:nt/2)'+dp); ...
         sin(2*pi*gamma*t(nt/2+1:end)')*0];
xa    = xa/A;
xb1   = xb1/A;
xb2   = xb2/A;
xg    = xg/A;

clear Psig
[~,~,~,Psig(:,1)] = preprocdata(xa , 'freq', f, 'w0', w0, 'dt', dt);
[~,~,~,Psig(:,2)] = preprocdata(xb1, 'freq', f, 'w0', w0, 'dt', dt);
[~,~,~,Psig(:,3)] = preprocdata(xb2, 'freq', f, 'w0', w0, 'dt', dt);
[~,~,~,Psig(:,4)] = preprocdata(xg , 'freq', f, 'w0', w0, 'dt', dt);
loglog(f,Psig,'k--');

%% Set Legend, Axes, and Title
if strcmp(meth,'PCC')
    h = legend('total (sine+noise)', 'incoherent', 'coherent', ...
        'vol.cond.', 'sine w/o noise');
end
% IC or PLI
if any(strcmp(meth,{'wPLI', 'PLI', 'IC'}))
    h = legend('total (sine+noise)', 'incoherent', 'coherent', ...
        'sine w/o noise');
end

set(h, 'Interpreter', 'Latex');
ylim([1e-8 1e-6])
set(gca, 'Ysca', 'log')
xlim([1 70])
xlabel('f [Hz]', 'Interpreter', 'Latex')
ylabel('PSD [V$^2$/Hz]', 'Interpreter', 'Latex')
title([meth ', $\omega_0=' num2str(w0) ', n_\sigma=' num2str(nsig) '$'], ...
    'Interpreter', 'Latex');

% Save and close
% print(fig1, fname, '-depsc')
% close(fig1)
set(fig1, 'Visi', 'On')


%% Reconstruct time series for PCC only at the moment
if strcmp(meth,'PCC')
    % Figure
    fig1 = figure('Papersize', [10 4], 'PaperPosition', [0.75 0.5 8.5 3], ...
            'PaperPositionmode', 'manual', 'Visible', 'on');
    
    % Define W_coh, W_inc and W_vc
    Ph  = angle(Wxy(:,:,1,2))/pi*180;
    coh = C(:,:,1,2)>sig & abs(Ph)>phase_thresh;
    inc = C(:,:,1,2)<=sig;
    vc  = C(:,:,1,2)>sig & abs(Ph)<=phase_thresh;

    % Define the modified matrices of wavelet coefficients as given in Equations (8)-(10)
    Wcoh = W1(:,1:ds:end,1);
    Wcoh(~coh) = 0;
    Winc = W1(:,1:ds:end,1);
    Winc(~inc) = 0;
    Wvc  = W1(:,1:ds:end,1);
    Wvc(~vc) = 0;

    % wavelet reconstruction
    [ts_coh, P_coh] = wave_recstr ( Wcoh, f, w0 );
    [ts_inc, P_inc] = wave_recstr ( Winc, f, w0 );
    [ts_vc,  P_vc]  = wave_recstr ( Wvc, f, w0 );
    [ts_tot, P_tot, varWT] = wave_recstr ( W1(:,1:ds:end,1), f, w0 );
    fprintf('Variance of original and reconstructed time series : %f/%f.\n', ...
        sigma2, varWT)
    fprintf('Error of %f%%.\n', (sigma2-varWT)/sigma2)

    plot(t, ts_vc+0.015, 'g-', t, ts_inc-0.015, 'r-', t, ts_coh, 'b-')
    xlim([0 T])
    ylim([-0.025 0.025]);
    xlabel('time [s]', 'interp', 'latex')
    ylabel('PCC signals [V]', 'interp', 'latex') 

    %% Plot original and reconstructed time series
    %% this only looks good if frequency resolution is sufficiently high
    % figure, plot(t,X(:,1), t,ts_tot)
end