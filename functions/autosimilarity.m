function [similarity, sems] = autosimilarity(clusters)
%AUTOSIMILARITY Compute autosimilarity of cluster sequence
% Marcus A. Triplett. (2018).

[num_trials, num_clusters] = size(clusters);
fprintf('Computing autosimilarity function for %i trials and %i clusters...\n', num_trials, num_clusters);
similarity = zeros(num_clusters - 1, 1);
sems = zeros(num_clusters - 1, 1);
counts = zeros(num_clusters - 1, 1);
for t = 1:num_trials
    for offset = 1:num_clusters - 1
        fprintf('calculating ASI for arg %i...\n', offset);
        for cc = 1:num_clusters - offset
            bm = best_match(clusters{t, cc}, clusters{t, cc + offset});
            similarity(offset) = similarity(offset) + bm;
            sems(offset) = sems(offset) + bm^2;
            counts(offset) = counts(offset) + 1;
        end
    end
    fprintf('Trial %i complete.\n', t);
end

% single pass mean and SEM calculation
for offset = 1:num_clusters - 1
    var = (sems(offset) - (similarity(offset)^2)/counts(offset))/(counts(offset) - 1);
    sems(offset) = sqrt(var/counts(offset));
    similarity(offset) = similarity(offset)/counts(offset);
end

similarity = [1; similarity];
sems = [0; sems];

