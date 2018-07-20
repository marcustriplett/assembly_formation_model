function sm = best_match(cluster_a, cluster_b)
% BEST_MATCH    Returns the best match similarity between two given clusterings

len_a = max(cluster_a);
len_b = max(cluster_b);

sets_a = cell(len_a, 1);
sets_b = cell(len_b, 1);
for i = 1:len_a
    sets_a{i} = find(cluster_a == i);
end
for j = 1:len_b
    sets_b{j} = find(cluster_b == j);
end

vec_dist = @(u, v) length(intersect(u, v))/length(union(u, v));

sm = 0;
for i = 1:len_a
    mx_a = 0; % d ranges from 0 to 1
    for j = 1:len_b
        % min over j: d(s_i, s_j)
        ds = vec_dist(sets_a{i}, sets_b{j});
        mx_a = max(ds, mx_a);
    end
    sm = sm + mx_a;
end

for j = 1:len_b
    mx_b = 0;
    for i = 1:len_a
        ds = vec_dist(sets_b{j}, sets_a{i});
        mx_b = max(ds, mx_b);
    end
    sm = sm + mx_b;
end

sm = sm/(len_a + len_b);