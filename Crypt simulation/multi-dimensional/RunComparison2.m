%% --- Configuration ---
ROWS  = 3;           % Number of layers
COLS  = 4;           % Cells per ring (Small enough for Exact Graph)
ALPHA = 10.0;        % Differentiation hierarchy strength
Mode = 'spatial';    % Simulation mode

% Lambda range (Relative Fitness Cost)
lambda_vals = 0.4 : 0.2 : 5.0; 

% Monte Carlo Settings
NUM_TRIALS = 3000; % Higher trials for better precision

% Result Storage
% Structure: results(lambda_idx, layer_idx, method_idx)
% Method 1: Exact Spatial (Graph Solver)
% Method 2: Mean Field Approximation (TheoreticalSolver Count State)
% Method 3: Monte Carlo Simulation
results = zeros(length(lambda_vals), ROWS, 3);

fprintf('--- Starting 3-Way Validation (Rows=%d, Cols=%d, Alpha=%.1f) ---\n', ROWS, COLS, ALPHA);

%% --- Main Loop ---
for i = 1:length(lambda_vals)
    lam = lambda_vals(i);
    fprintf('Processing Lambda = %.2f ... ', lam);
    
    % --- Method 1: Exact Spatial Graph Solver ---
    [probs_spatial, ~] = ExactGraphSolver(ROWS, COLS, lam, ALPHA);
    results(i, :, 1) = probs_spatial;
    
    % --- Method 2: Mean-Field Approximation Solver ---
    [probs_meanfield, ~] = calculateExactFixation(ROWS, COLS, lam, ALPHA);
    results(i, :, 2) = probs_meanfield;
    
    % --- Method 3: Monte Carlo Simulation ---
    for lay = 1:ROWS
        wins = 0;
        start_pos = [lay, 1]; % Fixed column 1
        
        parfor t = 1:NUM_TRIALS
            % no visualization
            [winID, ~] = HexCryptSimulation(start_pos, ROWS, COLS, lam, ALPHA, false, Mode);
            if winID == -1 
                wins = wins + 1;
            end
        end
        results(i, lay, 3) = wins / NUM_TRIALS;
    end
    
    fprintf('Done.\n');
end

%% --- Visualization 1: Probability Curves ---
fig1 = figure('Color', 'w', 'Position', [100, 400, 1200, 400]);
N_pop = ROWS * COLS;

% Calculate Standard Moran (Well-Mixed, ignoring alpha)
moran_curve = (1 - 1./lambda_vals) ./ (1 - (1./lambda_vals).^N_pop);
moran_curve(lambda_vals == 1) = 1 / N_pop;

for lay = 1:ROWS
    subplot(1, ROWS, lay);
    hold on;
    
    % 1. Exact Spatial (The Truth for Def 2.1)
    plot(lambda_vals, results(:, lay, 1), 'b-o', 'LineWidth', 2, 'MarkerSize', 4, ...
        'DisplayName', 'Exact Spatial (Matrix)');
    
    % 2. Simulation (Should overlay Blue)
    plot(lambda_vals, results(:, lay, 3), 'rx', 'LineWidth', 1.5, 'MarkerSize', 8, ...
        'DisplayName', 'Simulation (Monte Carlo)');
        
    % 3. Mean Field (Approximation)
    plot(lambda_vals, results(:, lay, 2), 'g--', 'LineWidth', 2, ...
        'DisplayName', 'Mean-Field (Approximation)');
    
    % 4. Standard Moran Reference
    plot(lambda_vals, moran_curve, 'k:', 'LineWidth', 1, ...
        'DisplayName', 'Standard Moran (No structure)');
    
    % Styling
    xlabel('Fitness Advantage (\lambda)');
    ylabel('Fixation Probability');
    title(sprintf('Start in Layer %d (Weight \\alpha^%d)', lay, ROWS-lay+1)); % Check alpha logic
    grid on;
    ylim([-0.05 1.05]); 
    
    if lay == 1
        legend('Location', 'northwest', 'FontSize', 8);
    end
    
    hold off;
end
sgtitle(sprintf('Validation: Fixation Probabilities (R=%d, C=%d, \\alpha=%.0f)', ROWS, COLS, ALPHA));


%% --- Visualization 2: Residual Analysis (Error vs Simulation) ---
% Residual = Theory - Simulation
% A value near 0 means perfect prediction.

fig2 = figure('Color', 'w', 'Position', [100, 100, 1200, 300]);

% Calculate absolute Errors
err_spatial = results(:, :, 1) - results(:, :, 3); % Exact - Sim
err_meanfield = results(:, :, 2) - results(:, :, 3); % MF - Sim

% Find max error range for consistent Y-axis
max_err = max([abs(err_spatial(:)); abs(err_meanfield(:))]);
ylim_range = [-max_err*1.1, max_err*1.1];

for lay = 1:ROWS
    subplot(1, ROWS, lay);
    hold on;
    
    % Zero Line
    yline(0, 'k-', 'LineWidth', 1);
    
    % Plot Spatial Residuals
    plot(lambda_vals, err_spatial(:, lay), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'DisplayName', 'Exact Matrix Residual');
    
    % Plot Mean Field Residuals
    plot(lambda_vals, err_meanfield(:, lay), 'g--s', 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'DisplayName', 'Mean-Field Residual');
        
    % Styling
    xlabel('Fitness Advantage (\lambda)');
    ylabel('Prediction Error (Theory - Sim)');
    title(sprintf('Residuals for Layer %d', lay));
    grid on;
    ylim(ylim_range);
    
    if lay == 1
        legend('Location', 'best', 'FontSize', 8);
    end
    
    hold off;
end

sgtitle('Residual Analysis: How much does Different Theory deviate from Simulation?');
drawnow;

%% --- Summary Statistics ---
% Calculate Sum of Squared Errors (SSE) across all trials and layers
sse_spatial = sum(err_spatial(:).^2);
sse_meanfield = sum(err_meanfield(:).^2);

% Calculate Mean Absolute Error (MAE)
mae_spatial = mean(abs(err_spatial(:)));
mae_meanfield = mean(abs(err_meanfield(:)));

fprintf('\n--- Performance Summary ---\n');
fprintf('Total Squared Error (Exact Matrix): %.6f  (MAE: %.6f)\n', sse_spatial, mae_spatial);
fprintf('Total Squared Error (Mean Field):   %.6f  (MAE: %.6f)\n', sse_meanfield, mae_meanfield);

if sse_spatial < sse_meanfield
    fprintf('>> Exact Matrix Solver is %.1fx more accurate than Mean-Field.\n', sse_meanfield/sse_spatial);
else
    fprintf('>> Mean-Field is unexpectedly more accurate (check simulation parameters).\n');
end