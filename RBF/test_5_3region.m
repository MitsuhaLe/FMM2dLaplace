flag_geom = 'star';
nContNodes = 100;
nCorners = 5;
C = contourData(nContNodes, flag_geom, nCorners);
X = zeros(nContNodes, 3);
Y = X;
for i = 1:3
    r = 1-i*0.25;
    X(:,i) = r*C(:,1);
    Y(:,i) = r*C(:,2);
end
interPts = [X(:), Y(:)];
nIntPts = 3 * nContNodes;

bdyPts = contourData(nContNodes, flag_geom, nCorners);
plot(bdyPts(:,1),bdyPts(:,2),'.r', ...
    interPts(:,1),interPts(:,2),'.k' )
legend('边界点','内部点')