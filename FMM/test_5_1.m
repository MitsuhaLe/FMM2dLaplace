%%% 5.1节
flag_geom = 'star';
nCorners = 5;
curvelen = 2*pi;
iprec = 5;

%%% plot the boundary and target nodes
n = 400;
C = contourData(n,flag_geom,nCorners);

nTargets = 1000;
rng("default");
rng(7)
ttint    = 2*pi*rand(1,nTargets);
rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
r = 0.7*rmin*rand(1,nTargets);
target    = [r .* cos(ttint); r .* sin(ttint)];
uu_ref = bdyfunc1(target);
figure(1)
plot(C(1,[1:n,1]), C(4,[1:n,1]    ),'k.-',...
     target(1,:), target(2,:),'b.');
legend('Contour C','Target points xxs')
axis equal

%%% 计算误差，分析收敛性
nNodes = 200:200:4000;
len = length(nNodes);
errmax = zeros(2,len);

FLAG = ['s2';'s6'];

for j = 1:2
    flag_pot = FLAG(j,:);
    for i = 1:len
        nContNodes = nNodes(i);
        disp(i)
        %%% Set the geometry.
        C = contourData(nContNodes,flag_geom,nCorners);
    
        %%% Compute the Dirichlet boundary condition.
        uu_dir = bdyfunc1(C([1,4],:));
    
        %%% Create the "correction" matrix.
        A_corr = fmmCorrect(C,flag_pot,curvelen);
    
        %%% Solve a reference problem.
        gmresTol = 1e-13;
        [sigma,~,~,gmIter] = gmres(@(x) forgmres(x, C, A_corr, flag_pot,iprec),uu_dir,20,gmresTol,40);
        uu = evalPot(target,C,curvelen,sigma,flag_pot );
        errmax(j,i) = max(abs(uu - uu_ref)); 
    end
end
figure(2)
semilogy(nNodes,errmax(1,:),'s-',nNodes,errmax(2,:),'d-')
legend('s2','s6')
xlabel('N')
ylabel(errmax)

%%% 计算矩阵向量积的时间
nNodes = 200:200:4000;
len = length(nNodes);
averageTime = zeros(2,len);
averageTimeDirect = zeros(2,len);

for j = 1:2
    flag_pot = FLAG(j,:);
    for i = 1:len
        nContNodes = nNodes(i);
        disp(i)
        %%% Set the geometry.
        C = contourData(nContNodes,flag_geom,nCorners);

        %%% Compute the Dirichlet boundary condition.
        uu_dir = bdyfunc1(C([1,4],:));

        %%% Create the "correction" matrix.
        A_corr = fmmCorrect(C,flag_pot,curvelen);

        %%% Solve a reference problem.
        gmresTol = 1e-13;
        tic
        [~,~,~,gmIter] = gmres(@(x) forgmres(x, C, A_corr, flag_pot,iprec),uu_dir,20,gmresTol,40);
        averageTime(j,i) = toc / (gmIter(1) * 20 + gmIter(2));
        tic
        [~,~,~,gmIter] = gmres(@(x) forgmres_direct(x, C, A_corr, flag_pot),uu_dir,20,gmresTol,40);
        averageTimeDirect(j,i) = toc / (gmIter(1) * 20 + gmIter(2));
    end
end

figure(3)
plot(nNodes,averageTime(1,:),'.-',nNodes,averageTime(2,:),'x-', ...
    nNodes,averageTimeDirect(1,:),'s-',nNodes,averageTimeDirect(2,:),'d-')
legend('s2','s6','s2direct','s6direct')
xlabel('N')
ylabel('Time (s)')