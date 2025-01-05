function x = forgmres_direct(x,C,A_corr,flag_pot)
%UNTITLED4 此处提供此函数的摘要
%   此处提供详细说明

nContNodes = size(C,2);
source = [C(1,:); C(4,:)];
target = source;
nSources = nContNodes;
nTargets = nContNodes;
h = (2*pi)/nContNodes;

sigma = A_corr * x;

% call rfmm 
% rfmm2dpart(iprec,nsource,source,ifcharge,
% charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,
% ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg)
% [U]=RFMM2DPART(IPREC,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC);

% initialization
charge = zeros(1,nContNodes);
dipvec = zeros(2,nContNodes);
dipstr = zeros(1,nContNodes);
speed = sqrt(C(2,:).*C(2,:) + C(5,:).*C(5,:));

if ((flag_pot == 's2') | (flag_pot == 's6') | (flag_pot == 'su')) % Laplace: single 
    ifcharge = 1;
    ifdipole = 0;
    charge = - h/(2*pi) * x' .* speed;
elseif (flag_pot == 'dr')
    ifcharge = 0;
    ifdipole = 1;
    nn1   =  C(5,:)./speed;
    nn2   = -C(2,:)./speed;
    dipstr = -h/(2*pi) *x' .* speed;
    dipvec = [nn1; nn2];
else
    disp("flag_pot is wrong")
end



% charge = zeros(1,nContNodes);

%  charge = ones(1,nContNodes);

ifpot = 1;
ifgrad = 0;
ifhess = 0;
ifpottarg = 0;
ifgradtarg = 0;
ifhesstarg = 0;


% tic
% iprec=4;
% r2dpartdirect
% [U]=rfmm2dpart(iprec,nSources,source,ifcharge,charge,ifdipole,dipstr, ...
%     dipvec,ifpot,ifgrad,ifhess,nTargets,target,ifpottarg,ifgradtarg,ifhesstarg);
[U]=r2dpartdirect(nSources,source,ifcharge,charge,ifdipole,dipstr, ...
    dipvec,ifpot,ifgrad,ifhess,nTargets,target,ifpottarg,ifgradtarg,ifhesstarg);
% total_time=toc
% speed=(nsource+ntarget)/total_time
x = reshape(U.pot,nContNodes,1) + sigma;


end