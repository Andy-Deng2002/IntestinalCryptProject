%% --- Configuration ---
ROWS  = 3;           % Number of layers
COLS  = 4;           % Cells per ring (Small enough for Exact Graph)
ALPHA = 10.0;        % Differentiation hierarchy strength
Mode = 'spatial';    % Simulation mode

% Lambda range (Relative Fitness Cost)
lambda_vals = 0.4 : 0.3 : 5.0; 

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
    
    % --- Method 2: Macroscopic Approximation Solver ---
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
% Layout Configuration
margin_l = 0.08;  % Left margin
margin_r = 0.02;  % Right margin
margin_b = 0.15;  % Bottom margin for x-label
margin_t = 0.12;  % Top margin for title
gap      = 0.00;  % No gap between plots

% Calculate geometry
total_w_avail = 1 - margin_l - margin_r;
sp_width = total_w_avail / ROWS;
sp_height = 1 - margin_b - margin_t;

% Adjust figure size to ensure subplots are roughly square
base_height = 500;
desired_aspect = 1; % Width / Height of one plot
% (sp_width * FigW) / (sp_height * FigH) = 1
% FigW = FigH * (sp_height / sp_width)
fig_width = base_height * (sp_height / sp_width) * 0.9; % Reduced by 15%

fig1 = figure('Color', 'w', 'Position', [100, 400, fig_width, base_height]);
N_pop = ROWS * COLS;

% Calculate Standard Moran (Well-Mixed, ignoring alpha)
moran_curve = (1 - 1./lambda_vals) ./ (1 - (1./lambda_vals).^N_pop);
moran_curve(lambda_vals == 1) = 1 / N_pop;

% Define colors for academic presentation
color_exact = [0 0.4470 0.7410];      % Standard Blue
color_mf    = [0.4660 0.6740 0.1880]; % Standard Green
color_sim   = [0.8500 0.3250 0.0980]; % Standard Orange

for lay = 1:ROWS
    % Calculate exact position [left bottom width height]
    pos_left = margin_l + (lay-1) * (sp_width + gap);
    
    subplot('Position', [pos_left, margin_b, sp_width, sp_height]);
    hold on;
    
    % 1. Exact Spatial (The Truth for Def 2.1)
    % Solid line with filled markers
    h1 = plot(lambda_vals, results(:, lay, 1), '-o', ...
        'Color', color_exact, 'LineWidth', 1.5, ...
        'MarkerSize', 5, 'MarkerFaceColor', color_exact, ...
        'DisplayName', 'MicSMP (Exact)');
    
    % 2. Mean Field (Approximation)
    % Dashed line
    h2 = plot(lambda_vals, results(:, lay, 2), '--', ...
        'Color', color_mf, 'LineWidth', 2, ...
        'DisplayName', 'MacSMM (Approximation)');
        
    % 3. Simulation (Monte Carlo)
    % Cross markers, no line (discrete data points)
    h3 = plot(lambda_vals, results(:, lay, 3), 'x', ...
        'Color', color_sim, 'LineWidth', 1.5, 'MarkerSize', 7, ...
        'DisplayName', 'Simulation (Monte Carlo)');
    
    % Calculate MAEs requested for text display
    mae_mf_exact_lay = mean(abs(results(:, lay, 2) - results(:, lay, 1)));
    mae_sim_exact_lay = mean(abs(results(:, lay, 3) - results(:, lay, 1)));
    
    % Display text on the graph
    txt_str = {sprintf('MAE (MacSMM): %.5f', mae_mf_exact_lay), ...
               sprintf('MAE (Monte Carlo): %.5f', mae_sim_exact_lay)};
               
    % Position text at Top Right, but slightly left to avoid border
    % Or Top Left is usually safer for fixation curves which start at 0
    text(0.05, 0.95, txt_str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1 1 1 0.8], ...
        'EdgeColor', 'none', ...
        'FontSize', 10, 'FontName', 'Arial');

    % Styling
    xlabel('Proliferative advantage (\lambda)', 'FontSize', 10, 'FontWeight', 'bold');

    % Adaptive Y-Axis: Only show label and ticks for the first graph
    if lay == 1
        ylabel('Fixation Probability', 'FontSize', 12, 'FontWeight', 'bold');
    else
        set(gca, 'YTickLabel', []); % Remove Y-axis tick numbers
    end
    
    title(sprintf('Start Layer %d', lay), 'FontSize', 12, 'FontWeight', 'bold');
    
    grid on;
    set(gca, 'Box', 'on', 'LineWidth', 1.2, 'FontSize', 11, 'FontName', 'Arial', ...
        'GridAlpha', 0.15);
    ylim([-0.05 1.05]); 
    xlim([min(lambda_vals)-0.1, max(lambda_vals)+0.1]);
    
    if lay == 1
        legend([h1, h2, h3], 'Location', 'southeast', 'FontSize', 9);
    end
    
    hold off;
end
sgtitle(sprintf('Fixation Probability Curves: (M=%d, N=%d, \\alpha=%.1f)', ROWS, COLS, ALPHA), ...
    'FontSize', 15, 'FontWeight', 'bold');


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
        'DisplayName', 'Exact Graph Residual');
    
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
fprintf('Total Squared Error (Exact Graph): %.6f  (MAE: %.6f)\n', sse_spatial, mae_spatial);
fprintf('Total Squared Error (Mean Field):   %.6f  (MAE: %.6f)\n', sse_meanfield, mae_meanfield);

if sse_spatial < sse_meanfield
    fprintf('>> Exact Matrix Solver is %.1fx more accurate than Mean-Field.\n', sse_meanfield/sse_spatial);
else
    fprintf('>> Mean-Field is unexpectedly more accurate (check simulation parameters).\n');
end