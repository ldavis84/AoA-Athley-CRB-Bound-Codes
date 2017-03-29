function [p,a] = jones2pol(u)
% [p,a] = jones2pol(u);
% Takes a complex 2-vector (Jones) and maps it onto the polarization
% parameters and overall gain a:
% u = a*[exp(1i*p(2)) * cos(p(1)); sin(p(1))];
% u       -- 2 x m set of vectors to convert into polarization params
% p       -- 2 x m real polarization parameters, each in (-180,180]
% a       -- 1 x m complex scalars needed to get original vectors

z = exp(1i*angle(u(2,:)));
v = ([1;1]*conj(z)) .* u;   

a = z .* sqrt(sum(abs(u).^2,1));

p2 = angle(v(1,:));
v(1,:) = v(1,:) .* exp(-1i*p2);
v = real(v);

p1 = atan2(v(2,:),v(1,:));

p = 180/pi*[p1; p2];

end
