function A_corr = fmmCorrect(C,flag_pot,curvelen)
%%% Extract various parameters.
nContNodes     = size(C,2);
h        = curvelen/nContNodes;  % step length

if strlength(flag_pot) ~= 2
    disp("flag_pot is wrong");
end

% Several different kernels are tested:
%   's2'   Laplace single layer potential with  2nd order K-R correction
%   's6'   Laplace single layer potential with  6th order K-R correction
%   'su'   Laplace single layer potential with 10th order K-R correction
%   'dr'   Laplace double layer potential + point source (smooth kernel) 
%%% 搞笑的是，这里他都选取两个字母表示，是为了后面条件判断简单。
%%% 我想改成 slpkr02 slpkr06 slpkr10 dlp反而不行，不过也可以再dlp这里加一个后缀，还是尊重原作者的符号吧。

%%% Check what order corrections are required,
%%% and prepare the relevant index vectors for A_corr.
if flag_pot == 's2'
    % （ii,jj） represent the row and column of diagonal entries
    MU2 = [ 0.7518812338640025 + 0.1073866830872157e1; ...
        -0.7225370982867850 - 0.6032109664493744];
    JJ = ones(4,1)*(1:nContNodes);
    II = [(-2):(-1),1:2]' * ones(1,nContNodes) + JJ;
    II = mod(II-1,nContNodes)+1;
    KR = MU2([2,1,1,2]) * ones(1,nContNodes);
    ii = II(:)';
    jj = JJ(:)';
    kr = KR(:)';
elseif flag_pot == 's6'
    MU6 = [ 0.2051970990601250e1 + 0.2915391987686505e1;...
        -0.7407035584542865e1 - 0.8797979464048396e1;...
        0.1219590847580216e2 + 0.1365562914252423e2;...
        -0.1064623987147282e2 - 0.1157975479644601e2;...
        0.4799117710681772e1 + 0.5130987287355766e1;...
        -0.8837770983721025   - 0.9342187797694916];
    JJ = ones(12,1)*(1:nContNodes);
    II = [(-6):(-1),1:6]' * ones(1,nContNodes) + JJ;
    II = mod(II-1,nContNodes)+1;
    KR = [MU6(end:(-1):1);MU6] * ones(1,nContNodes);
    ii = II(:)';
    jj = JJ(:)';
    kr = KR(:)';
elseif flag_pot == 'su'
    MU10 = [ 0.3256353919777872D+01 + 0.4576078100790908D+01;...
        -0.2096116396850468D+02 - 0.2469045273524281D+02;...
        0.6872858265408605D+02 + 0.7648830198138171D+02;...
        -0.1393153744796911D+03 - 0.1508194558089468D+03;...
        0.1874446431742073D+03 + 0.1996415730837827D+03;...
        -0.1715855846429547D+03 - 0.1807965537141134D+03;...
        0.1061953812152787D+03 + 0.1110467735366555D+03;...
        -0.4269031893958787D+02 - 0.4438764193424203D+02;...
        0.1009036069527147D+02 + 0.1044548196545488D+02;...
        -0.1066655310499552D+01 - 0.1100328792904271D+01];
    JJ = ones(20,1)*(1:nContNodes);
    II = [(-10):(-1),1:10]' * ones(1,nContNodes) + JJ;
    II = mod(II-1,nContNodes)+1;
    KR = [MU10(end:(-1):1);MU10] * ones(1,nContNodes);
    ii = II(:)';
    jj = JJ(:)';
    kr = KR(:)';
elseif flag_pot == 'dr'
else
    disp("flag_pot is wrong");
end

% jj长度超过C的列，但是下面的式子是指标索引
% nn是单位外法向量
% dd是求积节点的坐标


if ((flag_pot == 's2') | (flag_pot == 's6') | (flag_pot == 'su')) % Laplace: single 
    speed = sqrt(C(2,jj).*C(2,jj) + C(5,jj).*C(5,jj));  % |G'(t)|
    dd1   = C(1,jj) - C(1,ii);
    dd2   = C(4,jj) - C(4,ii);
    aa   = -(h/(4*pi))*log(dd1.*dd1 + dd2.*dd2).*speed;
    aakr   = kr.*aa;
elseif (flag_pot == 'dr')
    dg1  = C(2,:);
    dg2  = C(5,:);
    ddg1 = C(3,:);
    ddg2 = C(6,:);
    aakr = (h/(4*pi))*(-ddg2.*dg1 + ddg1.*dg2)./(dg1.^2 + dg2.^2) - ...
                0.5*ones(1,nContNodes);
else
    disp("flag_pot is wrong")
end


%指定非零元位置和元素，以及整个矩阵的大小
if flag_pot == 'dr'
    ii = 1:nContNodes;
    jj = ii;
end
A_corr = sparse(ii,jj,aakr,nContNodes,nContNodes);

end