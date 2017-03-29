function g = hpolgain(i,Az,El,lam)
% g = hpolgain(i,Az,El,lam)]
% 3 element array picks up the electric field due to the horizontal pol
% d = [cos(El)*cos(Az); cos(El)*sin(Az); sin(El)];
% u = [-sin(Az); cos(Az); 0];
% v = [-sin(El)*cos(Az); -sin(El)*sin(Az); cos(El)];

g = 0*Az;

i = rem(i-1,3) + 1;

if (i == 1)
    g = -sin(Az);
elseif (i == 2)
    g = cos(Az);
end
g = g *4;

end