function vv = evalPot(target, C, curvelen, sigma, flag_pot)
% After getting sigma, evaluate the function value on target Nodes.
%
% The contour Nodes in C are actually sources.



nContNodes  = size(C,2);
nTargets  = size(target,2);
h  = curvelen/nContNodes;

%%% Create index vectors.
ii = reshape((1:nTargets)' * ones(1,nContNodes),1,nTargets*nContNodes);
jj = reshape(ones(nTargets,1) * (1:nContNodes) ,1,nTargets*nContNodes);

%%% Compute some vectors 
speed = sqrt(C(2,jj).*C(2,jj) + C(5,jj).*C(5,jj));
nn1   =  C(5,jj)./speed;
nn2   = -C(2,jj)./speed;
dd1   =  C(1,jj) - target(1,ii);
dd2   =  C(4,jj) - target(2,ii);

if (flag_pot == 'dr')
  ddsq  = dd1.*dd1 + dd2.*dd2;
  ee    = (-h/(2*pi))*((nn1.*dd1 + nn2.*dd2)./ddsq).*speed; 
elseif ((flag_pot == 's2') | (flag_pot == 's6') | (flag_pot == 'su'))
  ee    = -(h/(4*pi))*log(dd1.*dd1 + dd2.*dd2).*speed;
end

EVAL = full(sparse(ii,jj,ee));
vv   = EVAL*sigma;

end