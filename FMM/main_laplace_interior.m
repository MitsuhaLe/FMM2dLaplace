function [gmIter,errmax,errmsq,nContNodes] = main_laplace_interior(cParams,bdyfunc,tol,iprec,solver,target)
%MAIN_LAPLACE_INTERIOR 2d interior laplace equation with dirichlet 
% boundary condtion. It's an integral equation based solver which uses 
% FMM and  singular quadrature -- Kapur-Rokhlin quadrature.
% 
% 
% Input parameters:
% 
% cParams parameters of boundary/Contour. Now only 'star' shape
% has been implemented.
% cParams = [nCorners, rmin] 
%       rmin  minimum distance between center and the boundary.
% 
%
% bdyfunc: accept function handle of exact solution 
%
% tol: precision, which decides the number of coutour nodes.
%  
% iprec - FMM precision flag
%
%     -2 => tolerance =.5d0   =>  
%     -1 => tolerance =.5d-1  =>  1 digit 
%      0 => tolerance =.5d-2  =>  2 digits
%      1 => tolerance =.5d-3  =>  3 digits
%      2 => tolerance =.5d-6  =>  6 digits
%      3 => tolerance =.5d-9  =>  9 digits
%      4 => tolerance =.5d-12 => 12 digits
%      5 => tolerance =.5d-15 => 15 digits
% 
% solver:
%   's2'   Laplace single layer potential with  2nd order K-R correction
%   's6'   Laplace single layer potential with  6th order K-R correction
%   'su'   Laplace single layer potential with 10th order K-R correction
%   'dr'   Laplace double layer potential + point source (smooth kernel) 
%
% target: target meshgrid, in other words, points that require evaluation
%


nCorners = cParams(1);
rmin = cParams(2);

%%% check whether there are input 'slover' and 'target'.
if nargin == 4
    solver = 'dr';
    nTargets = 1000;
    rng("default");
    rng(7)
    ttint    = 2*pi*rand(1,nTargets);
    r = 0.7*rmin*rand(1,nTargets);
    target    = [r .* cos(ttint); r .* sin(ttint)];
elseif nargin == 5
    nTargets = 1000;
    rng("default");
    rng(7)
    ttint    = 2*pi*rand(1,nTargets);
    r = 0.7*rmin*rand(1,nTargets);
    target    = [r .* cos(ttint); r .* sin(ttint)];    
end

%%% initialization
flag_pot = solver;
maxNodes = 1e4;
nContNodes = 200;
flag_geom = 'star';
curvelen = 2*pi;
errmax = 1;

%%% Compute the exact solution.
% bdyfunc = @bdyfunc;
uu_ref = bdyfunc(target);

% if bdyfunc == 'bdyfunc1'
%     nSources = 3;
%     charge = nrand(1,nSources);
%     ttext = 2*pi*rand(1,nSources);
%     rmax = sqrt(max(C(1,:).^2 + C(4,:).^2));
%     source = 1.6*[rmax*cos(ttext);rmax*sin(ttext)];
%     bdyfunc = @bdyfunc;
%     uu_dir = bdyfunc(C([1,4],:),source,charge);
%     uu_ref = bdyfunc(target);
% else
%     bdyfunc = @bdyfunc;
%     uu_dir = bdyfunc(C([1,4],:));  % C([1,4],:) 
%     uu_ref = bdyfunc(target);
% end

while errmax > tol && nContNodes <= maxNodes

    %%% Set the geometry.
    C = contourData(nContNodes,flag_geom,nCorners);

    %%% Compute the Dirichlet boundary condition.
    uu_dir = bdyfunc(C([1,4],:));

    %%% Create the "correction" matrix.
    A_corr = fmmCorrect(C,flag_pot,curvelen);

    %%% Solve a reference problem.
    gmresTol = 0.1*tol;
    [sigma,~,~,gmIter] = gmres(@(x) forgmres(x, C, A_corr, flag_pot,iprec),uu_dir,20,gmresTol,20);

    %%% Check the error 
    uu = evalPot(target,C,curvelen,sigma,flag_pot );
    errmax = max(abs(uu - uu_ref));
    errmsq = sqrt(sum((abs(uu - uu_ref).^2))/length(uu));
    
    nContNodes = nContNodes + 1000;
end
% disp(nContNodes);
% fprintf(1,'      %2s    %10.2e   %10.2e \n',flag_pot,errmax,errmsq)

%%% plot the boundary and target nodes
n = 400;
C = contourData(n,flag_geom,nCorners);
plot(C(1,[1:n,1]), C(4,[1:n,1]    ),'k.-',...
     target(1,:), target(2,:),'b.');
legend('Contour C','Target points xxs')
axis equal

end