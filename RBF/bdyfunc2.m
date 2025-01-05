function bdy = bdyfunc2(target)
% Exact solution with form exp(x) * sin(y)
% 

bdy = exp(target(1,:)) .* sin(target(2,:));
bdy = bdy';
end