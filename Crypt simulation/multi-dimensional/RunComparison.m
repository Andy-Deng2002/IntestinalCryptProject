%% --- Configuration ---
ROWS  = 2;           % Number of layers
COLS  = 16;          % Cells per ring
ALPHA = 10.0;        % Differentiation hierarchy strength
Mode = 'spatial';    % Simulation mode

% Lambda range (Relative Fitness Cost)
lambda_vals = 1.5 : 0.4 : 3.5; 

% Monte Carlo Settings
NUM_TRIALS = 1000;

% Result Storage
% Structure: results(lambda_idx, layer_idx, type_idx)
% Dimension 3: 1 = Theoretical, 2 = Simulation
results = zeros(length(lambda_vals), ROWS, 2);

fprintf('--- Starting Comparison (Rows=%d, Alpha=%.1f) ---\n', ROWS, ALPHA);

%% --- Main Loop ---
for i = 1:length(lambda_vals)
    lam = lambda_vals(i);
    fprintf('Processing Lambda = %.2f ... ', lam);
    
    % --- 1. Theoretical Exact Solution ---
    [probs_exact, ~] = calculateExactFixation(ROWS, COLS, lam, ALPHA);
    
    % Store exact results for all layers
    for lay = 1:ROWS
        results(i, lay, 1) = probs_exact(lay);
    end
    
    fprintf('[Theory Done] ');
    
    % --- 2. Monte Carlo Simulation ---
    % Loop through each starting layer dynamically
    for lay = 1:ROWS
        wins = 0;
        start_pos = [lay, round(COLS/2)]; % Start in middle of current layer
        
        % Use parfor for speed
        parfor t = 1:NUM_TRIALS
            [winID, ~] = HexCryptSimulation(start_pos, ROWS, COLS, lam, ALPHA, false, Mode);
            if winID == -1 
                wins = wins + 1;
            end
        end
        
        results(i, lay, 2) = wins / NUM_TRIALS;
    end
    
    fprintf('-> Sim Done.\n');
end

%% --- Visualization ---
figure('Color', 'w', 'Position', [100, 100, 400 * ROWS, 500]);

for lay = 1:ROWS
    subplot(1, ROWS, lay);
    hold on;
    
    % Plot Theoretical Curve
    plot(lambda_vals, results(:, lay, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'Theory');
    
    % Plot Simulation Points
    plot(lambda_vals, results(:, lay, 2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Simulation');
    
    % Styling
    xlabel('Adhesion (\lambda)');
    ylabel('Fixation Probability');
    
    title(sprintf('Start in Layer %d', lay));
    grid on;
    legend('Location', 'best');
    ylim([-0.05 1.05]); % Keep probability view tidy
    
    hold off;
end

sgtitle(sprintf('Fixation Probability: Theory vs Simulation (#Row=%d, \\alpha=%.0f, mode=%s)', ROWS, ALPHA, Mode));

drawnow;