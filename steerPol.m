function [b,gradb] = steerPol(egain,lam,Rarray,Euler,TrueAzEl,TruePol)
% [b,gradb] = steerPol(egain,lam,Rarray,Euler,TrueAzEl,TruePol);
% Finds the true normalized steering vector and its derivatives with
% respect to azimuth and elevation, and the two polarization parameters.
% Az is clockwise from north (bearing), El up from horizontal
%    El rotates up from xy plane to z, Az about z right-handed (CCW).
% Source direction is d = [cosd(El)*cosd(Az); cosd(El)*sind(Az); sind(El)];
% egain     -- function returning a single element's complex gain for one 
%              of two orthogonal polarizations.  If g() is
%              the function, one passes @g in the call to this function and
%              G1 = g(i,Azvals,Elvals,lam,ipol) must return the two
%              amplitude gains of element i in matrix G with same dims
%              as Azvals, Elvals, when these are in degrees.
%              ipol = 1 or 2 to select one or the other polarization.
% lam       -- wavelength of the energy
% Rarray    -- 3 x nelem locations of the array elements.   
% Euler     -- 3 x nelem Euler angles in deg indicating the orientation of
%              each element in the global system. First rotation is about
%              +z (Azimuth), then about rotated +x (Elev), then
%              about rotated-rotated +y (Roll). If empty, then all
%              rotations are assumed zero.  
% TrueAzEl  -- 1 x 2 true source azimuth, elevation in degrees.  
% TruePol   -- 1 x 2 true polarization Jones parameters p1 and p2 in deg
% b         -- nelem x 1 unit manifold vector for incoming wave including 
%              its polarization
% gradb     -- nelem x 4 gradient of b with respect to az, el, p1, p2
%              all units deg

nelem = size(Rarray,2);

Az = TrueAzEl(1);
El = TrueAzEl(2);

[mu,dmudp1,dmudp2] = pol2jones(TruePol);     % Correct Jones vector

% d = look direction, u local horizontal, v local vertical

normal = false;

[G1,G2] = arrayManifoldPol(egain,lam,Rarray,Euler,Az,El,normal);
G = [G1 G2];
b = G*mu;

bnorm = norm(b);

if (bnorm == 0)
    gradb = zeros(nelem,4);
    return;
end

b = b / bnorm;       % normalized

if (nargout == 1)
    return; 
end

delta = 1e-5;    % small enough angle perturbation for any array, deg

[G1pa,G2pa] = arrayManifoldPol(egain,lam,Rarray,Euler,Az+delta/2,El,normal);
[G1ma,G2ma] = arrayManifoldPol(egain,lam,Rarray,Euler,Az-delta/2,El,normal);

dG1dAz = (G1pa - G1ma)/delta;
dG2dAz = (G2pa - G2ma)/delta;

[G1pe,G2pe] = arrayManifoldPol(egain,lam,Rarray,Euler,Az,El+delta/2,normal);
[G1me,G2me] = arrayManifoldPol(egain,lam,Rarray,Euler,Az,El-delta/2,normal);

dG1dEl = (G1pe - G1me)/delta;
dG2dEl = (G2pe - G2me)/delta;

dbdAz_raw = [dG1dAz dG2dAz]*mu / bnorm;
dbdEl_raw = [dG1dEl dG2dEl]*mu / bnorm;

dbdp1_raw = G*dmudp1 / bnorm;
dbdp2_raw = G*dmudp2 / bnorm;

dbdAz = dbdAz_raw - b*real(b'*dbdAz_raw);
dbdEl = dbdEl_raw - b*real(b'*dbdEl_raw);

dbdp1 = dbdp1_raw - b*real(b'*dbdp1_raw);
dbdp2 = dbdp2_raw - b*real(b'*dbdp2_raw);

gradb = [dbdAz dbdEl dbdp1 dbdp2];

end