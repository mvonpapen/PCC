function [Pcoh, Pinc, Pvc] = pcc ( f, W, C, coi, sig, Wxy, phase_thresh )
%% Phase-coherence classification and subsequent estimation of avg PSD for each signal class
% 
% 
%   INPUT:
%           f:      frequency vector
%           W:      Wavelet coefficient matrix [Nf x Nt x Nch]
%           C:      Coherence matrix (can also be imaginary part of coherency)
%                   [Nf x Nt x Nch x Nch]
%           coi:    Cone of influence vector
%           sig:    Significance threshold according to nsig and w0
%           Wxy:    Cross-Wavelet coefficient matrix [Nf x Nt x Nch x Nch]
%           phase_thres: Phase threshold in degrees
% 
%   OUTPUT:
%           Pcoh    PSD of coherent signal [Nf x Nch x Nch]
%           Pinc    PSD of incoherent signal [Nf x Nch x Nch]
%           Pvc     PSD of volume-conducted signal [Nf x Nch x Nch]
%                   (=coherent signal with |dPhi|<phase_thresh)
% 
% Author: Michael von Papen
% Date: 22.04.16
%
% Notice: Please acknowledge the use of this program in any publications:
%   "Software for application of PCC was provided by von Papen et al. (2017)
%    and is available on github."
%
% Reference: von Papen, M. et al. (2017) ‘Phase-coherence classification: 
%            A new wavelet-based method to separate local field potentials 
%            into local (in)coherent and volume-conducted components’, 
%            Journal of Neuroscience Methods, 291, pp. 198–212. 
%            doi: 10.1016/j.jneumeth.2017.08.021.
%
% We would be grateful to receive a copy of such publications (m.von.papen@fz-juelich.de).


% Set phase threshold to determine volume-conducted signal
% Phi_c=15.5 deg according to w0=12 and nsigma=6
if nargin<8
    phase_thresh = 15.5;
end


if nargin<6 && nargout>2
    error('Cannot compute Pvc without Phase matrix!')
end


% Set significance to 1% threshold
if nargin<5
    sig = 0.41;
end


% Check sizes of input matrices
[nf, nt, nx]        = size ( W );
[nf2, nt2, nxw, ny] = size ( C );
if nf~=nf2 || nt~=nt2 || nx~=nxw || nx~=ny
    error('Sizes of W and C do not match!')
end
clear nf2 nt2 nxw ny


% Set cone of influence
if nargin<3
    coi = ones(nt,nx)*9e9;
end


% Pre-allocate variables
Pcoh = NaN(nf,nx,nx);
Pinc = NaN(nf,nx,nx);
if nargout>2
    Pvc = NaN(nf,nx,nx);
end
if nargin<6
    Ph = ones(size(C))+phase_thresh;
else
    Ph = abs(angle(Wxy)/pi*180);
end


%% Begin loop over channels
for i=1:nx
    
    Xw = W(:,:,i);
    Xc = coi(:,i)';
    
    for j=1:nx
        
        if i==j
            continue
        end
        
        % find indices for (in)coherent coefficients
        coh = C(:,:,i,j)>sig & Ph(:,:,i,j)>=phase_thresh ...
            & Ph(:,:,i,j)<=180-phase_thresh;
        inc = C(:,:,i,j)<=sig;
               
        % calculate PSD for coherent coefficients
        Wi = zeros(nf,nt);
        tmp = Xw;
        Wi(coh) = tmp(coh);
        Pcoh(:,i,j) = psdw( Wi, min([Xc; coi(:,j)']), f );
        
        % calculate PSD for incoherent coefficients
        Wi = zeros(nf,nt);
        tmp = Xw;
        Wi(inc) = tmp(inc);
        Pinc(:,i,j) = psdw( Wi, min([Xc; coi(:,j)']), f );
        
        % Optional: calculate PSD for coherent coefficients with dPhi ~ 0 & 180
        % corresponding to volume conduction
        if nargout>2
            vc = ~coh & ~inc;
            Wi = zeros(nf,nt);
            tmp = Xw;
            Wi(vc) = tmp(vc);
            Pvc(:,i,j) = psdw( Wi, min([Xc; coi(:,j)']), f ); 
        end          
        
    end
    
end     