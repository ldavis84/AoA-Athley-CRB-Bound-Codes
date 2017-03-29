function [psi,theta,phi] = Q2euler(Qgr)
% Qgr = euler2Q(psi,theta,phi);
% Converts direction cosine matrix to Euler angles.  Angles in degrees.
% First rotation is about +z (Azimuth), then about rotated -y (Elev), then
% about rotated-rotated +x (Roll).
% Qgr           -- 3x3 direction cosine matris, if vr is in the rotated
%                  system, then vg = Qgr * vr is in the original, and
%                  vr = Qgr' * vg.
% psi           -- azimuth, deg
% theta         -- elevation/pitch, deg
% phi           -- roll deg

theta = atan2d(Qgr(3,1), sqrt(Qgr(3,2)^2 + Qgr(3,3)^2));
sth = Qgr(3,1);

if (abs(1-abs(sth)) <= 100*eps)
    phi = 0;
    psi = atan2d(-Qgr(1,2),-Qgr(2,2));
else
    phi = atan2d(Qgr(3,2),Qgr(3,3));
    psi = atan2d(Qgr(2,1),Qgr(1,1));
end

end
