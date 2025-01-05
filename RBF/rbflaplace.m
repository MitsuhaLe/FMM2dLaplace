function [idx, errmax, errmsq] = rbflaplace(bdyfunc,nCorners)
%RBFLAPALCE MAIN_LAPLACE_INTERIOR 2d interior laplace equation 
% with dirichlet boundary condtion. 
%
% It's an collocation method with radius basis function 
% We refer to Gregory E. Fasshauer, page 358, program 39.2
% Kansa collocation for 2D elliptic eqn.
%
% The points in this code follow the form n x 2;
%



Rbf  = @(e,r) 1./sqrt(1+(e*r).^2);                                              % IMQ RBF
dxRbf = @(e,r,dx) -dx*e^2 ./ (1+(e*r).^2).^(3/2);                               % x方向一阶偏导数
dyRbf = @(e,r,dy) -dy*e^2 ./ (1+(e*r).^2).^(3/2);                               % y方向一阶偏导数
dxxRbf = @(e,r,dx) e^2 * ( 3*(e*dx).^2 -1-(e*r).^2 ) ./ (1+(e*r).^2).^(5/2);    % x方向二阶偏导数
dyyRbf = @(e,r,dy) e^2 * ( 3*(e*dy).^2 -1-(e*r).^2 ) ./ (1+(e*r).^2).^(5/2);    % x方向二阶偏导数


% Rbf  = @(e,r) exp(-(e*r).^2);                                                     % Gaussians                                   
% dxRbf = @(e,r,dx) exp(-(e*r).^2) .* (-2 * e^2 .* dx);                               % x方向一阶偏导数
% dyRbf = @(e,r,dy) exp(-(e*r).^2) .* (-2 * e^2 .* dy);                               % y方向一阶偏导数
% dxxRbf = @(e,r,dx) exp(-(e*r).^2) .* (-2 * e^2 .* dx).^2 -2 * e^2 * exp(-(e*r).^2);    % x方向二阶偏导数
% dyyRbf = @(e,r,dy) exp(-(e*r).^2) .* (-2 * e^2 .* dy).^2 -2 * e^2 * exp(-(e*r).^2);    % x方向二阶偏导数


shape = 0.5; % shapepara: the shape paramter e above

% initialization
flag_geom = 'star';
nContNodes = 50;
maxNodes = 400;
% nCorners = 5;
C = contourData(nContNodes, flag_geom, nCorners);
rmin = sqrt(min(C(:,1).^2 + C(:,2).^2));
nTargets = 100;
rng("default");
rng(7)
ttint = 2*pi*rand(nTargets,1);
r = 0.7*rmin*rand(nTargets,1);
target = [r .* cos(ttint), r .* sin(ttint)];

% Exact solution and its RHS
uu_ref = bdyfunc(target');

len = (maxNodes - nContNodes)/10 + 1;
errmax = zeros(1,len);
errmsq = zeros(1,len);
idx = zeros(1,len);

ilen = 1;
while nContNodes <= maxNodes

    % (equally spaced) boundary collocation points
    bdyPts = contourData(nContNodes, flag_geom, nCorners);
    
    % interior data sites (collocation points)
    C = bdyPts;
    X = zeros(nContNodes, 3);
    Y = X;
    for i = 1:3
        r = 1-i*0.25;
        X(:,i) = r*C(:,1);
        Y(:,i) = r*C(:,2);
    end
    interPts = [X(:), Y(:)];
    nIntPts = 3 * nContNodes;
    
    % ---------- Centers -------------------------------------------
    Ctrs = [ interPts; bdyPts ];  
    % ------------------------------------------------------------
    
    % laplace condition matrix \nable u = 0
    DM_intdata = DistanceMatrix( interPts, Ctrs );
    dx_int = DifferenceMatrix( interPts(:,1), Ctrs(:,1) );
    dy_int = DifferenceMatrix( interPts(:,2), Ctrs(:,2) );
    LCM = dxxRbf( shape, DM_intdata, dx_int )...
        + dyyRbf( shape, DM_intdata, dy_int );
    
    % boundary condition matrix u|_gamma = f(x)
    DM_bdydata = DistanceMatrix( bdyPts, Ctrs );
    BCM = Rbf( shape, DM_bdydata );
    
    CM  = [ LCM; BCM ];
    uu_dir = bdyfunc(bdyPts');
    rhs = [ zeros(nIntPts,1); uu_dir ];
%     [U,s,V] = csvd(CM);
%     [reg_corner1,~,~,~] = l_curve(U,s,rhs,'tsvd');
%     [ Coeff_tsvd, ~, ~ ] = tsvd(U,s,V,rhs,reg_corner1) ;
    Coeff = CM \ rhs;
%     Coeff = Coeff_tsvd;
    % evaluation and check errors
    DM_eval = DistanceMatrix( target, Ctrs );
    EM = Rbf( shape, DM_eval );
    uu = EM * Coeff;
    errmax(ilen) = max(abs(uu - uu_ref));
    errmsq(ilen) = sqrt(sum((abs(uu - uu_ref).^2))/length(uu));
    idx(ilen) = nContNodes;
%     fprintf('MAX error:  %e\n',errmax);
%     fprintf('RMS error:  %e\n',errmsq);

    nContNodes = nContNodes + 10;
    ilen = ilen + 1;
end




% figure, PF = zeros(Neval); PF(:) = Pf; surf(PF)
% figure, EXACT = zeros(Neval); EXACT(:) = exact; surf(EXACT)