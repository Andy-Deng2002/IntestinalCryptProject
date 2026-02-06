ROWS = 2;
COLS = 16;
ALPHA = 10.0;
LAMBDA = 2;

for i = 1:2
    % --- Target position ---
    % [Row, Col]
    % Row 1 = Base (Stem), Row N = Top (Lumen)
    TARGET_POS = [i, round(COLS/2)]; 
    
    % ---  ---statistical setting
    NUM_TRIALS = 5000;
    VISUALIZE = false;
    
    fprintf('Runnng Monte Carlo simulation for Fixation Probability...\n');
    fprintf('Grid: %dx%d | Alpha: %.2f | Lambda: %.2f\n', ROWS, COLS, ALPHA, LAMBDA);
    fprintf('Start Position: Row %d, Col %d\n', TARGET_POS(1), TARGET_POS(2));
    fprintf('Total Trials: %d\n', NUM_TRIALS);
    
    mutant_wins = 0;
    
    tic;
    
    parfor t = 1:NUM_TRIALS
    % TODO: learn parallel
        % run one simulation
        [winner, steps] = HexCryptSimulation(TARGET_POS, ROWS, COLS, LAMBDA, ALPHA, false, 'spatial');
        
        if winner == -1
            mutant_wins = mutant_wins + 1;
        end
    end
    
    elapsed_time = toc;
    fixation_prob = mutant_wins / NUM_TRIALS;
    
    % --- output ---
    fprintf('\n================ RESULTS ================\n');
    fprintf('Mutant Fixations: %d / %d\n', mutant_wins, NUM_TRIALS);
    fprintf('Fixation Probability: %.5f\n', fixation_prob);
    fprintf('Theoretical Neutral Prob (1/N): %.5f\n', 1/(ROWS*COLS));
    fprintf('Time Elapsed: %.2f seconds (%.2f ms/trial)\n', elapsed_time, (elapsed_time/NUM_TRIALS)*1000);
    fprintf('=========================================\n');
end