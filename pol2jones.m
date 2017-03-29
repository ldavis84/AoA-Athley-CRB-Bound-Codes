function [u,dudp1,dudp2] = pol2jones(p)
% [u,dudp1,dudp2] = pol2jones(p);
% Takes a real poloarization specification parameters and maps them 
% onto a complex 2-vector (Jones):
% p      -- 2 x m sets of real polarization specification parameters, p(1)
%           specifies fraction of each polarization, p(2) phase lead of 
%           first one.  angles in deg
% u      -- 2 x m Jones vectors 
%           = [exp(1i*pi/180*p(2)) * cosd(p(1)); sind(p(1))];
% dudp1  -- 2 x m derivatives of the vectors wrt first parameter
% dudp2  -- 2 x m derivatives of the vectors wrt 2nd parameter

if (length(p) == 2)
    p = p(:);
end

u = [exp(1i*pi/180*p(2,:)) .* cosd(p(1,:)); sind(p(1,:))];

if nargout == 1
    return;
end

m = size(p,2);

dudp1 = [-exp(1i*pi/180*p(2,:)) .* sind(p(1,:)); cosd(p(1,:))]*pi/180;
dudp2 = [1i*exp(1i*pi/180*p(2,:)) .* cosd(p(1,:)); zeros(1,m)]*pi/180;

end