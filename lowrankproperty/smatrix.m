function S = smatrix(x, y)
% single potential matrix
%   targets：x_i in R^2  (2, n), n>=2;
%   sources: y_j in R^2  (2, m), m>=2;
% well-separated condition. The purpose is to verify the low rank property.
% S_{ij} = ln | x_i - y_j |

xlen = length(x);
ylen = length(y);
S = zeros(xlen, ylen);
for ix = 1:xlen
    % temp = x(:, ix) - y;
    % temp = vecnorm(x(:, ix) - y);
    S(ix, :) = log(vecnorm(x(:, ix) - y));
end

end