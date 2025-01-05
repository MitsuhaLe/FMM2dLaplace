%%% 5.2节
% 两种问题区域，两种精确解，四种组合

%% Omega_1 + u_1
clc
clear

nContNodes = 200;
flag_geom = 'star';
nCorners = 5;
C = contourData(nContNodes,flag_geom,nCorners);
rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
cParams = [nCorners, rmin];

tol = 1e-13;
iprec = 5;

bdyfunc = @bdyfunc1;
FLAG = ['s2';'s6';'su';'dr'];
for i = 1:4
    solver = FLAG(i,:);
    [gmIter,errmax,errmsq,nContNodes] = main_laplace_interior(cParams,bdyfunc,tol,iprec,solver);
    fprintf(1,'      %2s  %2d %2d   %10.2e   %10.2e  %10d \n ', ...
        solver,gmIter(1),gmIter(2),errmax,errmsq,nContNodes-200);
end
disp(1)

%% Omega_1 + u_2
clear

nContNodes = 200;
flag_geom = 'star';
nCorners = 5;
C = contourData(nContNodes,flag_geom,nCorners);
rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
cParams = [nCorners, rmin];

tol = 1e-13;
iprec = 5;

bdyfunc = @bdyfunc2;
FLAG = ['s2';'s6';'su';'dr'];
for i = 1:4
    solver = FLAG(i,:);
    [gmIter,errmax,errmsq,nContNodes] = main_laplace_interior(cParams,bdyfunc,tol,iprec,solver);
    fprintf(1,'      %2s  %2d %2d   %10.2e   %10.2e  %10d \n ', ...
        solver,gmIter(1),gmIter(2),errmax,errmsq,nContNodes-200);
end
disp(2)

%% Omega_2 + u_1
clear

nContNodes = 200;
flag_geom = 'star';
nCorners = 3;
C = contourData(nContNodes,flag_geom,nCorners);
rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
cParams = [nCorners, rmin];

tol = 1e-13;
iprec = 5;

bdyfunc = @bdyfunc1;
FLAG = ['s2';'s6';'su';'dr'];
for i = 1:4
    solver = FLAG(i,:);
    [gmIter,errmax,errmsq,nContNodes] = main_laplace_interior(cParams,bdyfunc,tol,iprec,solver);
    fprintf(1,'      %2s  %2d %2d   %10.2e   %10.2e  %10d \n ', ...
        solver,gmIter(1),gmIter(2),errmax,errmsq,nContNodes-200);
end
disp(3)

%% Omega_2 + u_2
clear

nContNodes = 200;
flag_geom = 'star';
nCorners = 3;
C = contourData(nContNodes,flag_geom,nCorners);
rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
cParams = [nCorners, rmin];

tol = 1e-13;
iprec = 5;

bdyfunc = @bdyfunc2;
FLAG = ['s2';'s6';'su';'dr'];
for i = 1:4
    solver = FLAG(i,:);
    [gmIter,errmax,errmsq,nContNodes] = main_laplace_interior(cParams,bdyfunc,tol,iprec,solver);
    fprintf(1,'      %2s  %2d %2d   %10.2e   %10.2e  %10d \n ', ...
        solver,gmIter(1),gmIter(2),errmax,errmsq,nContNodes-200);
end
disp(4)
