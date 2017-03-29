function Qcrb = CRBAoA(egain,SNRdB,lam,Rarray,Euler,M,...
    TrueAzEl,TruePol)
% Qcrb = CRBAoA(egain,SNRdB,lam,Rarray,Euler,M,TrueAzEl,TruePol);
% Finds the CRB for an array of point apertures for a single signal
% when the covariance of the signal white Gaussian signal is unknown but
% the covariance of the noise is known.  Az CCW from x direction, El up
% from horizontal towards z.
%    El rotates up from xy plane to z, Az about z right-handed (CCW).
% Source direction is d = [cos(El)*cos(Az); cos(El)*sin(Az); sin(El)];
%    Errors are in terms of local du,dv where uhat is unit direction vector 
% along cross(z,d), where z = [0;0;1], i.e., a local azimuthal rotation.  
% dAz = du / cos(El).  dEl = dv (dv is small rotation up).  So 
% vhat = cross(d,uhat) is local elevation, and [d uhat vhat]
% form a right-handed local coordinate system.  The Jones vector is:
% pJ = [cosd(p1) .* exp(j*pi/180*p2); sind(p1)], where [p1; p2]
% are the Jones parameters in degrees
% egain     -- function returning a single element's complex gain for one 
%              of two orthogonal polarizations.  If g() is
%              the function, one passes @g in the call to this function and
%              G1 = g(i,Azvals,Elvals,lam,ipol) must return the two
%              amplitude gains of element i in matrix G with same dims
%              as Azvals, Elvals, when these are in degrees.
%              ipol = 1 or 2 to select one or the other polarization.
% SNRdB     -- 1 x nsnr the SNR values for analysis.
%              SNR is with respect to a the whole array 1 snapshot,
%              i.e., SNRpower = qx * norm(A)^2 / sig^2, where qx is
%              the signal m.s. at a gain-of-1 antenna, norm(A) is the norm
%              of the vector of gains of all apertures, and sig = noise rms
% lam       -- wavelength of the energy
% Rarray    -- 3 x nelem locations of the array elements.   
% Euler     -- 3 x nelem Euler angles in deg indicating the orientation of
%              each element in the global system. First rotation is about
%              +z (Azimuth), then about rotated +x (Elev), then
%              about rotated-rotated +y (Roll). If empty, then all
%              rotations are assumed zero.  
% M         -- number of sample points assumed
% TrueAzEl  -- 1 x 2 true source azimuth, elevation in degrees.  
% TruePol   -- 1 x 2 true polarization Jones parameters p1 and p2 in deg
% Qcrb      -- 4 x 4 x nsnr the minimal covariance matrix for the local
%              azimuth and elevation

[~,nelem] = size(Rarray);

nsnr = length(SNRdB);
Qmax = 180^2 * diag([1/3, 1/6]);

% Obtain the array unit pointing and pointing derivative.

[bo,gradbo] = steerPol(egain,lam,Rarray,Euler,TrueAzEl,TruePol);
gradbo = gradbo(:,1:2);

P = eye(nelem) - bo*bo';    % projection matrix

gPg = real(gradbo'*P*gradbo);

diagload = max(abs(gPg(:))) * 10*eps;
gPg = gPg + diagload*eye(2);

Fo = inv(gPg);   % Array strength factor, in terms of true az angle change
c = cosd(TrueAzEl(1));

Fo(1,:) = Fo(1,:)*c;   % Correct to local dAz
Fo(:,1) = Fo(:,1)*c;

% Now complete the CRB using the SNR, maximum is if you don't know anything
% Qmax = [1 0; 0 0.25] * pi^2/3;

S = 10.^(SNRdB/10);    % power SNR
Ratio = 0.5*(1+S) ./ S.^2 / M;

Qcrb = kron(Ratio,Fo);
Qcrb = reshape(Qcrb,2,2,nsnr);

for i = 1:nsnr
    
    % find minimum in definiteness sense
    
    Qtmp = mindef(Qcrb(:,:,i),Qmax);
    Qcrb(:,:,i) = Qtmp;
end

end