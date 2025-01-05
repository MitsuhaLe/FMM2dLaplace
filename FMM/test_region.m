%%% 两种问题区域边界
nContNodes = 1000;
flag_geom = 'star';
nCorners = 5;
C = contourData(nContNodes,flag_geom,nCorners);
figure(1)
plot(C(1,:),C(4,:),'k.')

nContNodes = 1000;
flag_geom = 'star';
nCorners = 3;
C = contourData(nContNodes,flag_geom,nCorners);
figure(2)
plot(C(1,:),C(4,:),'k.')

