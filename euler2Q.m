function Qgr = euler2Q(psi,theta,phi)
% Qgr = euler2Q(psi,theta,phi);
% Converts Euler angles to direction cosine matrix converting rotated
% coordinates to original.  Angles in degrees.
% First rotation is about +z (Azimuth), then about rotated -y (Elev), then
% about rotated-rotated +x (Roll).
% psi           -- azimuth, deg
% theta         -- elevation/pitch, deg
% phi           -- roll deg
% Qgr           -- 3x3 direction cosine matris, if vr is in the rotated
%                  system, then vg = Qgr * vr is in the original, and
%                  vr = Qgr' * vg.

cpsi = cosd(psi);
spsi = sind(psi);

cth = cosd(theta);
sth = sind(theta);

cphi = cosd(phi);
sphi = sind(phi);

Qgr = [cpsi*cth, -spsi*cphi - cpsi*sth*sphi, spsi*sphi - cpsi*sth*cphi; ...
    spsi*cth, cpsi*cphi - spsi*sth*sphi, -cpsi*sphi - spsi*sth*cphi; ...
    sth, cth*sphi, cth*cphi];

end