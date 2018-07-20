function M = set_diag_zero(M)

for i = 1:size(M, 1)
   M(i, i) = 0; 
end

