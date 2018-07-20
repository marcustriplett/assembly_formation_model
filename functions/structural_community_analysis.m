function [cluster, structural_data] = structural_community_analysis(M, structural_data, t)

cluster = cluster_jl(M, 0, 0, 0, 0); % third party community detection

structural_data.num_clusters(t) = length(cluster.SIZE{1, 1});
structural_data.avg_cluster_size(t) = mean(cluster.SIZE{1,1});
structural_data.modularity(t) = cluster.MOD;
structural_data.coeff_variation(t) = std(cluster.SIZE{1,1})/mean(cluster.SIZE{1,1});
