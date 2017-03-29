function [ainv,u,v] = pseudo(a,toler)
% [ainv,u,v] = pseudo(a,toler)
% finds the pseudo inverse of a via svd. u and v are the orthogonal bases
% for the column and row spaces of a, respectively.

[m,n] = size(a);
flip = n > m;

if (flip)
  [v,s,u] = svd(a',0);
  s = diag(s);

  if (s(1) == 0)
    ainv = zeros(n,m);
    return;
  end
  
  ii = find(s >= toler*s(1));
  s = s(ii).^(-1);
  v = v(:,ii);
  u = u(:,ii);
  
  ainv = v*((s(:)*ones(1,m)).*u');
else
  [u,s,v] = svd(a,0);
  s = diag(s);

  if (s(1) == 0)
    ainv = zeros(n,m);
    return;
  end
  
  ii = find(s >= toler*s(1));
  s = s(ii).^(-1);
  u = u(:,ii);
  v = v(:,ii);
  
  ainv = (v.*(ones(n,1)*s(:)'))*u';
end
