% ASSEMBLY_FORMATION_NETWORK Simulates a binary neuron model of assembly formation from
%   initially random synaptic connectivity.
%   
%   Marcus A. Triplett. (2018).
%
%   See corresponding article Emergence of spontaneous assembly activity
%   in developing neural networks without afferent input. 
%   M. A. Triplett et al. (2018). To appear.

addpath functions

tic
num_trials = 1;

% drawing
draw_on = 1;
pause_on = 0; % enable for regular pausing.
save_on = 0;

% set trial duration
duration = 120000;
pause_interval = 5000;

% viewing window size
window_size = 5000;

% clustering data
clustering_interval = 1000;
sample_cluster_time = 100000; % starting time for calculating autosimilarity properties
if not(mod(duration, clustering_interval) == 0)
    error('clustering_interval must divide duration');
end

% set population parameters
Nn = 125; % total population size
Ne = ceil(0.80 * Nn); % excitatory population size
Ni = Nn - Ne; %inhibitory population size

% activity parameters
mu = 0.002; % mean background firing rate
sdv = 0.0003; % standard dev of background firing rate
background_rate = max(normrnd(mu, sdv, [Nn, 1]), 0);
activation_threshold = 0.1;
threshold = activation_threshold * ones(Nn, 1);

% plasticity and stimulation
eta = 0.04; % learning rate

dur_ov_cl = duration/clustering_interval;
hour = 3600;
minute = 60;

% trial means
tm_structural_data = struct();
som_structural_data = struct();
structural_fields = {'num_clusters', 'avg_cluster_size', 'modularity', 'event_freqs', 'coeff_variation'};
structural_field_titles = {'Num clusters', 'Avg cluster size', 'Modularity', 'Event freq', 'CV'};
for fn = 1:numel(structural_fields)
    tm_structural_data.(char(structural_fields(fn))) = zeros(dur_ov_cl, 1);
    som_structural_data.(char(structural_fields(fn))) = zeros(dur_ov_cl, 1);
end

saved_clusters = cell(num_trials, floor((duration - sample_cluster_time)/clustering_interval));

% create directories to save simulation data
if save_on
    mkdir(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f', eta, mu, sdv, activation_threshold, num_trials), '0.', ''))
    mkdir(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/weight_matrices', eta, mu, sdv, activation_threshold, num_trials), '0.', ''))
    mkdir(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/summary_statistics', eta, mu, sdv, activation_threshold, num_trials), '0.', ''))
    mkdir(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/simulations', eta, mu, sdv, activation_threshold, num_trials), '0.', ''))
end

% begin trial loop
for trial = 1:num_trials

% reset params
x_e = zeros(Ne, 1);
x_i = zeros(Ni, 1);
activities = zeros(Nn, duration);
mu_e = zeros(Ne, 1);

% generate and normalise weight matrices
W_ee = set_diag_zero(rand(Ne));
W_ei = normalise_rows(rand(Ne, Ni));
W_ie = normalise_rows(rand(Ni, Ne));
W_ii = normalise_rows(set_diag_zero(rand(Ni, Ni)));

% presynaptic and postsynaptic excitatory normalisation
for ii = 1:size(W_ee, 1)
    W_ee(ii, ii) = 0;
end
W_ee = W_ee ./ sum(W_ee, 1);
W_ee = W_ee ./ sum(W_ee, 2);

cluster = cluster_jl(W_ee, 0, 0, 0, 0);

% within-trial data
for fn = 1:numel(structural_fields)
    structural_data.(char(structural_fields(fn))) = zeros(dur_ov_cl, 1);
end

% file handling
if save_on
    for tx = [string('W_ei'), string('W_ie'), string('W_ii')]
        save(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/weight_matrices/weights_%s_trial_%.f', ...
            eta, mu, sdv, activation_threshold, num_trials, tx, trial), '0.', ''), char(tx))
    end
end

position_flag = 1;

% main simulation loop
for t = 1:duration
    x_e_temp = heaviside_origin(W_ee*x_e - W_ei*x_i - threshold(1:Ne) + (1 + activation_threshold) * binornd(ones(Ne, 1), background_rate(1:Ne)));
    x_i = heaviside_origin(W_ie*x_e - W_ii*x_i - threshold(Ne + 1:Nn) + (1 + activation_threshold) * binornd(ones(Ni, 1), background_rate(Ne + 1:Nn)));
    x_e = x_e_temp;

    activities(:, t) = [x_e; x_i];

    % update and renormalise excitatory weights
    mu_e = 1/t * ((t - 1) * mu_e + x_e);
    W_ee = max(W_ee + eta * (x_e - mu_e) * (x_e - mu_e)', 0); % covariance plasticity rule
    
    for ii = 1:size(W_ee, 1)
        W_ee(ii, ii) = 0;
    end
    W_ee = W_ee ./ sum(W_ee, 1);
    W_ee = W_ee ./ sum(W_ee, 2);
    
    % plot intermediate activity
    if mod(t, pause_interval) == 0
        [cluster, structural_data] = structural_community_analysis(W_ee, structural_data, t/clustering_interval);        
        ordering = order_by_cluster(Ne, cluster);
        if draw_on
            activities_sorted = [activities(ordering, max(1, t - window_size):t); activities(Ne + 1:Nn, max(1, t - window_size):t)];

            figure(1); clf;
            hold on;
            xlim([0, window_size]);
            ylim([0, Ne + 1]);
            for i = 1:Ne
                ac = find(activities_sorted(i, :));
                scatter(ac, i*ones(1, length(ac)), 8, 'black', 'o',  'filled');
            end
            box on;
            set(gca, 'linewidth', 3, 'YTick', 0:10:Ne, 'fontsize', 20, 'xtick', 0:1000:window_size, 'xticklabel', t-window_size:1000:t)
            ylabel('Neuron'); xlabel('Time');           
            set(gcf, 'color', 'white')
            if position_flag set(gcf, 'position', [388 488 560 420]); end;

            figure(2); clf;
            ordering_flipped = flipud(ordering);
            imagesc(W_ee(ordering_flipped, ordering_flipped));
            box on;
            set(gca, 'linewidth', 3, 'XTick', [], 'YTick', [], 'fontsize', 20);
            set(gcf, 'color', 'white');
            xlabel('Neuron'); ylabel('Neuron');
            colorbar;
            axis square;
            if position_flag set(gcf, 'position', [950 488 560 420]); position_flag = 0; end;

            drawnow;
            if pause_on, pause(); end;
        end
    end
%     
    if mod(t, clustering_interval) == 0
        structural_data.event_freqs(t/clustering_interval) = mean(1/clustering_interval * sum(activities(1:Ne, max(1, t - clustering_interval):t), 2));
        
        % compute structural community properties
        [cluster, structural_data] = structural_community_analysis(W_ee, structural_data, t/clustering_interval);
        structural_data.coeff_variation(t/clustering_interval) = std(cluster.SIZE{1,1})/mean(cluster.SIZE{1,1});
        
        ordering = order_by_cluster(Ne, cluster);
        
        % sample cluster
        if t == sample_cluster_time
            sample_cluster = cluster;
        end
        
        % compute community similarity properties
        if t >= sample_cluster_time
            saved_clusters{trial, (t - sample_cluster_time)/clustering_interval + 1} = cluster.COM{1,1};
        end
        
        if save_on
            save(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/weight_matrices/weights_W_ee_trial_%.f_t_%.f', ...
                eta, mu, sdv, activation_threshold, num_trials, trial, t), '0.', ''), 'W_ee')
        end
        
        
    if mod(t, 5000) == 0
       fprintf('(Trial %d) %d steps completed of %d \n', trial, t, duration)
    end
    
    end
end

if save_on
    save(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/simulations/activities_trial_%.f', eta, mu, sdv, activation_threshold, num_trials, trial), '0.', ''), ...
        'activities', 'structural_data', 'saved_clusters', 'W_ee', 'W_ei', 'W_ie', 'W_ii')
end

% update mean and variance parameters
for j = 1:numel(structural_fields)
    f = structural_fields{j};
    tm_structural_data.(f) = tm_structural_data.(f) + structural_data.(f);
    som_structural_data.(f) = som_structural_data.(f) + structural_data.(f).^2;
end

end % end trial loop

toc

% divisor for unbiased variance estimator
if num_trials > 1
    var_div = num_trials - 1;
else
    var_div = 1;
end

% tm_s still hold sums, calculate variance
sem_structural_data = struct();
for fn = 1:numel(structural_fields)
    sem_structural_data.(char(structural_fields(fn))) = zeros(dur_ov_cl, 1);
end

for j = 1:numel(structural_fields)
    f = structural_fields{j};
    var_estimator = (som_structural_data.(f) - (tm_structural_data.(f).^2)/num_trials)/var_div;
    sem_structural_data.(f) = sqrt(var_estimator/num_trials);
end

% divide through by num_trials to obtain means
for j = 1:numel(structural_fields)
    f = structural_fields{j};
    tm_structural_data.(f) = tm_structural_data.(f)/num_trials;
end

% autosimilarity already accounts for num_trials
[tm_autosim, sem_autosim] = autosimilarity(saved_clusters);
tm_structural_data.('autosimilarity') = tm_autosim;
sem_structural_data.('autosimilarity') = sem_autosim;

% plotting axis
xaxis_cl = duration/dur_ov_cl:duration/dur_ov_cl:duration;

if save_on
    % save structures
    tm_structural_data.('x_axis') = xaxis_cl';
    sem_structural_data.('x_axis') = xaxis_cl';
    
    save(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/summary_statistics/structural_avgs', eta, mu, sdv, activation_threshold, num_trials), '0.', ''), 'tm_structural_data');
    save(strrep(sprintf('bnet_sim_eta_%0.2f_mu_%0.3f_sig_%0.4f_theta_%.2f_ntrial_%0.0f/summary_statistics/structural_sems', eta, mu, sdv, activation_threshold, num_trials), '0.', ''), 'sem_structural_data'); 
end

if draw_on
    figure(3);
    for j = 1:numel(structural_fields)
        subplot(3, 2, j);
        f = structural_fields{j};
        if num_trials > 1
            errorbar(tm_structural_data.(f), sem_structural_data.(f),  'linewidth', 2)
        else
            plot(tm_structural_data.(f), 'linewidth', 2);
        end
        set(gcf, 'color', 'white');
        set(gca, 'fontsize', 10, 'linewidth', 2);
        xlim([0, numel(tm_structural_data.(f))]);
        box off; axis square;
        xlabel('Sampling timepoint'); ylabel(structural_field_titles{j});
    end
    
    subplot(3, 2, numel(structural_fields) + 1);
    plot(tm_autosim, 'k-.', 'linewidth', 2);
    set(gcf, 'color', 'white')
    set(gca, 'fontsize', 10, 'linewidth', 2, 'xtick', [0, numel(tm_autosim)], 'xticklabel', [0, numel(tm_autosim)*clustering_interval])
    xlim([0, numel(tm_autosim)]);
    ylim([0, 1.05])
    xlabel('Time'); ylabel('Autosimilarity');
    box off; axis square;
end
