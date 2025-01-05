
% dsites: M-by-s
% ctrs: N-by-s
% DM(i,j) = || dsites_i - ctrs_j ||, M-by-N.

function DM = DistanceMatrix( dsites, ctrs )

[M,s] = size(dsites);
[N,s] = size(ctrs);
DM = zeros(M,N);
for d = 1 : s
    [ dr, cc ] = ndgrid( dsites(:,d), ctrs(:,d) );
    DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);