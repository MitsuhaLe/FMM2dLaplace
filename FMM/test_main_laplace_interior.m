%% bdyfunc1
clc

nContNodes = 200;
flag_geom = 'star';
nCorners = 5;
C = contourData(nContNodes,flag_geom,nCorners);
rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
cParams = [nCorners, rmin];

tol = 1e-10;
iprec = 5;
% solver = 'dr';

bdyfunc = @bdyfunc1;
FLAG = ['s2';'s6';'su';'dr'];
for i = 1:4
    solver = FLAG(i,:);
    [gmIter,errmax,errmsq,nContNodes] = main_laplace_interior(cParams,bdyfunc,tol,iprec,solver); 
    fprintf(1,'      %2s  %2d %2d   %10.2e   %10.2e  %10d \n ',solver,gmIter(1),gmIter(2),errmax,errmsq,nContNodes);
end

%% bdyfunc2 
clc

nContNodes = 200;
flag_geom = 'star';
nCorners = 5;
[C,~] = contourData(nContNodes,flag_geom,nCorners);
rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
cParams = [nCorners, rmin];

tol = 1e-10;
iprec = 5;
solver = 'dr';

bdyfunc = @bdyfunc2;
[errmax,errmsq] = main_laplace_interior(cParams,bdyfunc,tol,iprec,solver);

