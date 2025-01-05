function bdy = exactPot(target,source,charge)
% This function represents dirichlet boundary condition.
% 
%
% source:x_j  target:y_i charge:q_j 
% Compute sum_{j=1}^n -log||y_i - x_j||
%

nSources = length(charge);
nTargets = size(target,2);

% DD1_ij = y_i(1) - x_j(1)
% DD2_ij = y_i(2) - x_j(2)
DD1  = target(1,:)'*ones(1,nSources) - ones(nTargets,1)*source(1,:);   
DD2  = target(2,:)'*ones(1,nSources) - ones(nTargets,1)*source(2,:);

AA = -(1/(2*pi))*0.5*log(DD1.*DD1 + DD2.*DD2); 
bdy = AA*reshape(charge,nSources,1);

end
