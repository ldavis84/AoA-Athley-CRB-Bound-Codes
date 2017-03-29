function [Lo,gradL,hessL] = AoAderiv(egain,z,lam,Rarray,Euler,Azo,Elo,U)
% [Lo,gradL,hessL] = AoAderiv(egain,z,lam,Rarray,Euler,Azo,Elo,U);
% Computes likelihood measure L and its derivatives with respect to Az and
% El at a particular point using finite differences.
%    El rotates up from xy to z, Az about z right-handed (CCW).
% Source direction is d = [cos(El)*cos(Az); cos(El)*sin(Az); sin(El)];
%   The user provides a function G = elemgain(i,Azvals,Elvals,lam); This 
% function provides the complex gain of element i relative to the gain
% assumed in the SNR in a matrix conformal to Azvals, Elvals.  A single 
% polarization is  assumed.  
% egain     -- function returning a single element's complex gain for each
%              value of azimuth and elevation.  If g() is
%              the function, one passes @g in the call to this function and
%              G = g(i,Azvals,Elvals,lam,ipol) must return the complex
%              amplitude gain of element i in matrix G with same dimensions
%              as Azvals, Elvals, when these are in radians.
% z         -- 2 x 1 complex Jones vector giving the polarization in terms
%              of the vertical and horizontal polarizations, resp.
% lam       -- wavelength of the energy
% Rarray    -- 3 x nelem locations of the array elements
% Azo       -- Az value at which to evaluate, deg
% Elo       -- El value at which to evaluate, deg
% Euler     -- 3 x nelem Euler angles in deg indicating the orientation of
%              each element in the global system. First rotation is about
%              +z (Azimuth), then about rotated +x (Elev), then
%              about rotated-rotated +y (Roll). If empty, then all
%              rotations are assumed zero.  
% U         -- nelem x nelem matrix to form L = norm(U'*A(az,el))^2, 
%              where A(az,el) is unit antenna manifold vector at az,el
% Lo        -- value of the likelihood measure
% gradL     -- 1 x 2 gradient of L wrt Az, El in rad
% hessL     -- 2 x 2 Hessian matrix of L

delta = 1e-7*180/pi;   % angle increment to use for perturbations

Azvals = Azo + delta*[0 1 -1 -1 1 2 0 -2 0];
Elvals = Elo + delta*[0 1 1 -1 -1 0 2 0 -2];

% Evaluate L at points.

normalize = true;
A = arrayManifold(egain,z,lam,Rarray,Euler,Azvals,Elvals,normalize);
A = squeeze(A);

Lvals = sum(abs(U'*A).^2,1);

% Value

Lo = Lvals(1);

% Gradients

dLdAz = (Lvals(2) + Lvals(5) - Lvals(3) - Lvals(4)) / (4*delta);
dLdEl = (Lvals(2) + Lvals(3) - Lvals(4) - Lvals(5)) / (4*delta);
gradL = [dLdAz dLdEl];

% Hessian Matrix

d2LdAz2 = ((Lvals(6) - Lvals(1)) - (Lvals(1) - Lvals(8))) / (2*delta)^2;
d2LdEl2 = ((Lvals(7) - Lvals(1)) - (Lvals(1) - Lvals(9))) / (2*delta)^2;

d2LdAzdEl = ((Lvals(2) - Lvals(3)) - (Lvals(5) - Lvals(4))) / (2*delta)^2;

hessL = [d2LdAz2 d2LdAzdEl; d2LdAzdEl d2LdEl2];

end
