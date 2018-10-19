function [x, y] = synth_data( t, nl, A, dp )
%% SYNTH_DATA generates synthetic time series
%
% Synthetic time series created here are used as test for the PCC algorithm.
% The time series reflect some basic characteristics of LFP data.
% The time series show alpha, beta and gamma oscillations that vary their
% amplitude in time (switched on or off). Powerlaw noise with a 1/f decay
% and a SNR similar to typical LFP observations is added.

if nargin<4
    dp = 30; # delta phi
end
if nargin<3
    A = 1e3; # normalization factor
end
if nargin<2
    nl = 3; # noise level
end

%% Parameters
alpha = 10;
beta  = 20;
gamma = 50;
dt    = t(2)-t(1); # delta t, time resolution
nt    = length(t); # sample length

%% Synthetic time series
% Generate noise vectors
n1 = nl * powlawnoise(nt,1)';
n2 = nl * powlawnoise(nt,1)';
% Generate segments with alpha, beta and gamma oscillations
xa = [zeros(nt/5,1); sin(2*pi*alpha*(1:3*nt/5)'*dt); zeros(nt/5,1)];
xb = [sin(2*pi*beta*(1+(1:nt/2)'/nt*0.25).*(1:nt/2)'*dt); ...
      sin(2*pi*beta*(1+(nt/2+1:nt)'/nt*0.25).*(nt/2+1:nt)'*dt-dp)];
xg = [sin(2*pi*gamma*t(1:nt/2)'+dp); sin(2*pi*gamma*t(nt/2+1:end)')];
% Add noise and oscillations
x  = xa + xb + xg + n1;
% Repeat process for second time series
ya = [zeros(nt/5,1); sin(2*pi*alpha*(1:3*nt/5)'*dt+dp); zeros(nt/5,1)];
yb = sin(2*pi*beta*(1+t'/(nt*dt)*0.25).*t');
yg = [sin(2*pi*gamma*t(1:nt/2)'); sin(2*pi*gamma*t(nt/2+1:end)')];
y  = ya + yb + yg + n2;
% Normalize output to typical LFP magnitudes (mV)
x  = x / A;
y  = y / A;

end
