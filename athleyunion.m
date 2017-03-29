function [Q,Pi] = athleyunion(SNRdB,nsnap,SLLdB,SLdu,Qcrb)
% [Q,Pi] = athleyunion(SNRdB,nsnap,SLLdB,SLdu,Qcrb);
% Computes the Athley union bound for a stochastic Gaussian white signal
% on an array with given sidelobe levels and various SNRs.  Produces mean
% squared error estimates in dB of direction error (from boresight
% to true target).
% SNRdB    -- SNR value at which to compute the bound, 
%             SNR is with respect to a the whole array single snapshot,
%             i.e., SNRpower = qx * norm(A)^2 / sig^2, where qx is
%             the signal m.s. at a gain-of-1 antenna, norm(A) is the norm
%             of the vector of gains of all apertures, and sig = noise rms
% nsnap    -- number of snapshots or samples
% SLLdB    -- 1 x nsl levels of the sidelobes relative to the main beam, dB
%             should all be <= 0
% SLdu     -- dim x nsl direction error vector for each sidelobe
% Qcrb     -- ndim x ndim Cramer-Rao bound on MSE of the direction cosine
%             vector. 2 x 2 for planar case, scalar for linear, etc. 
% Q        -- ndim x ndim similar to Qcrb but includes threshold errors due
%             to sidelobes
% Pi       -- estimated probabilities for the sidelobes

NSL = length(SLLdB);    % number of sidelobes

NSR1 = 10^(-SNRdB / 10);      % noise to signal ratio, single snap

SLL = 10.^(SLLdB / 10);

[SLL,isort] = sort(SLL,'descend');      % sort biggest first
SLdu = SLdu(:,isort);

% Compute nominal conditional probability of each sidelobe 

% Ntmp = 2*nsnap-1;
NSRtmp = 4*NSR1 .* (1 + NSR1);

Pi = zeros(1,NSL);

warnid = 'MATLAB:nchoosek:LargeCoefficient';  % turn off warn nchoosek
warning('off',warnid);

for i = 1:NSL
    
    r2 = SLL(i);
    
    qtmp = sqrt(1 + NSRtmp / (1 - r2));
    q = (qtmp + 1) ./ (qtmp - 1);
    
   
    for m = 0:nsnap-1
        Pi(i) = Pi(i) + (q/(1+q)).^m * nchoosek(nsnap-1+m,m);
    end
    
    Pi(i) = Pi(i) ./ (1 + q).^nsnap;
end

warning('on',warnid);   % restore nchoosek warning

% Correct probabilities so that they are independent and all still < 1.  If
% they are small then same as Athley.  Assumption is, a lower sidelobe is
% considered ONLY if ALL higher sidelobes are not taken.

logtmp = [0 log(1 - Pi)];
OneminusP = exp(cumsum(logtmp));

Pi = Pi .* OneminusP(1:NSL);

% Compute Po and MSE contribution for each sidelobe at each SNR.
% Athley modified to include Qcrb also for all sidelobes--resulting in 
% Qcrb * 1 in the total.

Q = Qcrb;    % CRB obtains as a variance on all sidelobes too.

for i = 1:NSL
    
    du2 = SLdu(:,i)*SLdu(:,i)';
    Q = Q + Pi(i) * du2;
    
end

end