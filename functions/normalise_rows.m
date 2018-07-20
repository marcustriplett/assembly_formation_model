function M = normalise_rows(M)

for j = 1:size(M, 1)
   s = sum(M(j, :));
   if s ~= 0
       M(j, :) = 1/sum(M(j, :)) * M(j, :);
   end
end