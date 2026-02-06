% filepath: c:\Users\xiluo\Desktop\Crypt simulation\multi-dimensional\RunComparison3.m
%% --- Configuration ---
ROWS  = 4;           % Taller cylinder to see diffusion effects
COLS  = 8;           
ALPHA = 1.0;         % Milder alpha to allow more vertical movement

% We compare two modes using the same Monte Carlo engine
modes = {'spatial', 'well_mixed'};
colors = {'b-o', 'r--x'};
names = {'Spatial (MicSMP)', 'Well-Mixed (MacSMM)'};

% Parameters
lambda_vals = 0.5 : 0.5 : 5.0; 
NUM_TRIALS = 2000; 

% Storage: results(lambda_idx, mode_idx)
prob_results = zeros(length(lambda_vals), 2);
time_results = zeros(length(lambda_vals), 2);

fprintf('--- TIME DYNAMICS COMPARISON (Rows=%d, Cols=%d, Alpha=%.1f) ---\n', ROWS, COLS, ALPHA);
fprintf('Hypothesis: Probabilities match, but Spatial is SLOWER than Well-Mixed.\n\n');

%% --- Main Simulation Loop ---
for m = 1:2
    current_mode = modes{m};
    fprintf('Running Mode: %s ...\n', current_mode);
    
    for i = 1:length(lambda_vals)
        lam = lambda_vals(i);
        
        wins = 0;
        total_time_steps = 0;
        success_count = 0;
        
        % Always start at Bottom Layer (Stem Cell) for consistency
        start_pos = [ROWS, 1]; 
        
        % Parallel Loop for speed
        parfor t = 1:NUM_TRIALS
            % Run simulation (Visualize=false)
            [winID, steps] = HexCryptSimulation(start_pos, ROWS, COLS, lam, ALPHA, false, current_mode);
            
            if winID == -1 
                wins = wins + 1;
                % Only count time for SUCCESSFUL fixations (standard practice)
                % Conditioned on fixation occurring.
                total_time_steps = total_time_steps + steps;
                success_count = success_count + 1;
            end
        end
        
        % Store Probability
        prob_results(i, m) = wins / NUM_TRIALS;
        
        % Store Conditional Mean Fixation Time
        if success_count > 0
            time_results(i, m) = total_time_steps / success_count;
        else
            time_results(i, m) = NaN; % No fixations observed clearly
        end
        
        fprintf('  Lambda %.1f: Prob=%.2f, MeanSteps=%.0f\n', lam, prob_results(i,m), time_results(i,m));
    end
end

%% --- Visualization ---
figure('Color', 'w', 'Position', [100, 200, 1200, 450]);

% 1. Fixation Probability (The Check)
subplot(1, 2, 1); hold on; grid on;
for m = 1:2
    plot(lambda_vals, prob_results(:, m), colors{m}, 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', names{m});
end
xlabel('Proliferative Advantage (\lambda)');
ylabel('Fixation Probability');
title(sprintf('Fixation Probability'));
legend('Location', 'best');
ylim([-0.05 1.05]);

% 2. Time to Fixation (The Proof of Non-Lumpability)
subplot(1, 2, 2); hold on; grid on;
for m = 1:2
    plot(lambda_vals, time_results(:, m), colors{m}, 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', names{m});
end
xlabel('Proliferative Advantage (\lambda)');
ylabel('Mean Fixation Time (Steps)');
title(sprintf('Time to Fixation'));
legend('Location', 'best');

sgtitle(sprintf('Spatial vs. Well-Mixed Dynamics (M=%d, N=%d, \\alpha=%.1f)', ROWS, COLS, ALPHA));

fprintf('\nDONE. Check Figure 1.\n');
if mean(time_results(:,1)) > mean(time_results(:,2))
    fprintf('CONCLUSION: Spatial model is SLOWER. System is NOT Lumpable (Dynamics differ).\n');
else
    fprintf('CONCLUSION: Times are similar. System might be effectively Lumpable due to small size.\n');
end