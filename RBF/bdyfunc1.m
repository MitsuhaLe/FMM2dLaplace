function bdy = bdyfunc1(target)
% Exact solution with form sum_{j=1}^n - q_j log||y_i - x_j||
% source:x_j  charge:q_j 
% target or evaluation nodes:y_i
% 
% 

nSources = 3;
source = [2.5, -2, 0;
          0  ,  0, 2 ];
charge = ones(1,nSources);
nTargets = size(target,2);

%%% DD1_ij = y_i(1) - x_j(1), DD2_ij = y_i(2) - x_j(2)
DD1  = target(1,:)'*ones(1,nSources) - ones(nTargets,1)*source(1,:);   
DD2  = target(2,:)'*ones(1,nSources) - ones(nTargets,1)*source(2,:);

AA = -(1/(2*pi))*0.5*log(DD1.*DD1 + DD2.*DD2); 
bdy = AA*reshape(charge,nSources,1);

end