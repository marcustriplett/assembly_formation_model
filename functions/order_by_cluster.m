function ordering = order_by_cluster(N, cluster)
ordering = zeros(N, 1);
nc = length(cluster.SIZE{1, 1});
st = 1;
for i = 1:nc
    cl = find(cluster.COM{1, 1} == i)';
    ordering(st:st + length(cl) - 1) = cl;
    st = st + length(cl);
end