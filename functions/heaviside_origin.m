% Heaviside step function with origin set to zero.

function hout = heaviside_origin(x)
n = length(x);
hout = zeros(n, 1);
for i = 1:n
    if x(i) < 0
        hout(i) = 0;
    else
        hout(i) = 1;
    end
end