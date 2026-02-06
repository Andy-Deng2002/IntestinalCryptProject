% filepath: c:\Users\xiluo\Desktop\Crypt simulation\HexCryptSimulation.m
function [winningCloneID, timeToMonoclonal] = HexCryptSimulation(mutant_pos, rows, cols, mutant_lambda, alpha, visualize, simulation_mode)
% HexCryptSimulation: 
% Simulates clonal competition on a hexagonal grid.
%
% Inputs:
%   ...
%   simulation_mode: 'spatial' (default) or 'well_mixed' (matches theoretical mean-field)
%
% Interactive Controls: Space=Pause, Up=Speed Up, Down=Slow Down

    % -- Parameter Defaults --
    if nargin < 7, simulation_mode = 'spatial'; end % 'spatial' or 'well_mixed'
    if nargin < 6, visualize = true; end 
    if nargin < 5, alpha = 2; end
    if nargin < 4, mutant_lambda = 10; end
    if nargin < 3, cols = 16; end
    if nargin < 2, rows = 3; end 
    if nargin < 1, mutant_pos = [1, round(cols/2)]; end 

    N = rows * cols;
    
    % Initialize Grid
    % WT ID: 1 to N
    % Mutant ID: -1
    grid = reshape(1:N, rows, cols);
    grid(mutant_pos(1), mutant_pos(2)) = -1; 
    
    maxIter = 1e6; 
    time_step = 0;
    winningCloneID = NaN;

    % --- Visualization Initialization (Only if asked) ---
    hFig = []; hScatterWT = []; hScatterMut = [];
    X_plot = []; Y_plot = [];
    
    if visualize
        % Create Figure and bind key press events
        hFig = figure('Name', sprintf('Hex Crypt Sim [%s] (Space:Pause)', simulation_mode), ...
                      'Color', 'w', 'KeyPressFcn', @keyPressHandler, 'NumberTitle', 'off');
        ax = axes(hFig);
        hold(ax, 'on');
        
        % Precompute Hexagonal Coordinates (Hex Shift Logic)
        [R, C] = ndgrid(1:rows, 1:cols);
        X_plot = C;
        % Shift odd rows to the right by 0.5 (Spatial mode only visual preference, keep consistent)
        X_plot(mod(R, 2) == 1) = X_plot(mod(R, 2) == 1) + 0.5;
        Y_plot = R; 
        
        % Ensure Y-axis is normal (Low values at bottom, High values at top)
        set(ax, 'YDir', 'normal'); 
        
        % Fix Color Map Range
        colormap(ax, parula(N)); 
        clim(ax, [1 N]); 
        
        % Create persistent graphics objects
        hScatterWT = scatter(ax, [], [], 300, [], 'filled', 'Marker', 'o'); 
        hScatterMut = scatter(ax, [], [], 300, 'r', 'filled', 'Marker', 'h', 'MarkerEdgeColor', 'k');
        
        axis(ax, 'equal', 'off');
        title(ax, 'Initializing...');
        
        % Labels
        mx = max(X_plot(:));
        text(ax, mx + 0.5, 1, 'Base (Stem)', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
        text(ax, mx + 0.5, rows, 'Top (Lumen)', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
        
        % Interactive State
        guiState.paused = false;
        guiState.delay = 0.03; 
        set(hFig, 'UserData', guiState);
        
        % Initial Draw
        updateGraphics(hScatterWT, hScatterMut, grid, X_plot, Y_plot, time_step, guiState.delay);
    end

    while time_step < maxIter
        % ---------------------------------------------------------
        % CHECK WIN CONDITION
        % ---------------------------------------------------------
        has_mutant = any(grid(:) == -1);
        has_wt = any(grid(:) ~= -1);
        
        if ~has_mutant
            winningCloneID = 1; % WT Won
            break;
        elseif ~has_wt
            winningCloneID = -1; % Mutant Won
            break;
        end
        
        time_step = time_step + 1;
        
        % ---------------------------------------------------------
        % VISUALIZATION & INTERACTION
        % ---------------------------------------------------------
        if visualize
            if ~isvalid(hFig), visualize = false; continue; end
            guiState = get(hFig, 'UserData');
            
            while guiState.paused
                pause(0.1);
                if ~isvalid(hFig), visualize = false; break; end
                guiState = get(hFig, 'UserData');
                title(hFig.CurrentAxes, 'PAUSED');
            end
            if ~visualize, break; end 
            
            if guiState.delay > 0.005 || mod(time_step, 20) == 0
                updateGraphics(hScatterWT, hScatterMut, grid, X_plot, Y_plot, time_step, guiState.delay);
            end
            
            if guiState.delay > 0, pause(guiState.delay); end
        end
        
        % ---------------------------------------------------------
        % STEP 1: Proliferate (Select cell to proliferate based on adhesion)
        % ---------------------------------------------------------
        row_indices = (rows:-1:1)';
        layer_weights = alpha .^ row_indices;
        weight_matrix = repmat(layer_weights, 1, cols);
        
        % high lambda -> high chance to be selected to proliferate
        
        is_mutant = (grid == -1); 
        weight_matrix(is_mutant) = weight_matrix(is_mutant) * (mutant_lambda);
        
        all_weights = weight_matrix(:);
        % if sum(all_weights) == 0, all_weights = ones(size(all_weights)); end
        
        proliferate_idx_linear = randsample(rows*cols, 1, true, all_weights);
        [r_out, c_out] = ind2sub([rows, cols], proliferate_idx_linear);
        
        % ---------------------------------------------------------
        % STEP 2: REPLACEMENT (Selection of neighbor)
        % ---------------------------------------------------------
        
        if strcmp(simulation_mode, 'well_mixed')
            % --- WELL MIXED MODE (Mean Field) ---
            % Matches TheoreticalSolver logic: 
            % Probability of sampling from layer L is proportional to link strength (k_link)
            
            [~, k_same, k_up, k_down] = get_hex_connectivity(r_out, rows);
            
            % Neighbors are defined by layers, not specific cells
            weights = [];
            candidate_layers = [];
            
            % Same layer
            if k_same > 0
                weights(end+1) = k_same; 
                candidate_layers(end+1) = r_out;
            end
            % Up layer
            if k_up > 0
                weights(end+1) = k_up;
                candidate_layers(end+1) = r_out - 1;
            end
            % Down layer
            if k_down > 0
                weights(end+1) = k_down;
                candidate_layers(end+1) = r_out + 1;
            end
            
            % Pick a source layer
            chosen_layer_idx = randsample(length(candidate_layers), 1, true, weights);
            source_layer = candidate_layers(chosen_layer_idx);
            
            % Pick a random cell from that layer
            % CRITICAL: If source is same layer, cannot pick self!
            while true
                rand_col = randi(cols);
                if source_layer == r_out && rand_col == c_out
                    if cols > 1
                        continue; % Resample if picked self
                    else
                        % Edge case: only 1 column, cannot replace self from same layer
                        % Skip replacement or break logic? In theory, rate is 0.
                        % For now, just keep old value (do nothing)
                        rand_col = c_out; 
                        break; 
                    end
                else
                    break;
                end
            end
            
            grid(source_layer, rand_col) = grid(r_out, c_out);
            
        else
            % --- SPATIAL MODE (Standard Hex Grid) ---
            targets = getHexNeighbors(r_out, c_out, rows, cols);
            
            if ~isempty(targets)
                n_targets = size(targets, 1);
                chosen_idx = randsample(n_targets, 1);
                neighbor_pos = targets(chosen_idx, :);
                grid(neighbor_pos(1), neighbor_pos(2)) = grid(r_out, c_out);
            end
        end
    end
    
    if isnan(winningCloneID)
        % If loop finished maxIter without winner
         has_wt = any(grid(:) > 0);
         if ~has_wt, winningCloneID = -1; else, winningCloneID = 1; end
    end
    
    timeToMonoclonal = time_step;
    
    if visualize && isvalid(hFig)
       updateGraphics(hScatterWT, hScatterMut, grid, X_plot, Y_plot, time_step, 0);
       resStr = 'WT WON';
       if winningCloneID == -1, resStr = 'MUTANT WON'; end
       title(hFig.CurrentAxes, sprintf('DONE! %s', resStr));
    end
end

% ---------------------------------------------------------
% Helper: Graphics Update
% ---------------------------------------------------------
function updateGraphics(hWT, hMut, grid, X, Y, step, delay)
    mask_mut = (grid == -1); 
    mask_wt = ~mask_mut;
    
    if any(mask_wt(:))
        set(hWT, 'XData', X(mask_wt), 'YData', Y(mask_wt), 'CData', grid(mask_wt));
    else
        set(hWT, 'XData', [], 'YData', []);
    end
    
    if any(mask_mut(:))
        set(hMut, 'XData', X(mask_mut), 'YData', Y(mask_mut));
    else
        set(hMut, 'XData', [], 'YData', []);
    end
    
    title(get(hWT, 'Parent'), sprintf('Step: %d | Mutants: %d', step, sum(mask_mut(:))));
    drawnow limitrate;
end

% ---------------------------------------------------------
% Helper: Keyboard
% ---------------------------------------------------------
function keyPressHandler(src, event)
    state = get(src, 'UserData');
    switch event.Key
        case 'space', state.paused = ~state.paused; 
        case 'uparrow', state.delay = max(0, state.delay - 0.01);
        case 'downarrow', state.delay = state.delay + 0.05;
    end
    set(src, 'UserData', state);
end

% ---------------------------------------------------------
% Helper: Spatial Neighbors
% ---------------------------------------------------------
function neighbors = getHexNeighbors(r, c, maxR, maxC)
    left = c - 1;  if left < 1, left = maxC; end
    right = c + 1; if right > maxC, right = 1; end
    % Horizontal (Same layer)
    potential_neighbors = [r, left; r, right];
    
    % Odd row logic (shifted right) vs Even row logic
    if mod(r, 2) == 1 
        if r > 1 % Check Up
           c_up_right = c + 1; if c_up_right > maxC, c_up_right = 1; end
           potential_neighbors = [potential_neighbors; r-1, c; r-1, c_up_right];
        end
        if r < maxR % Check Down
           c_down_right = c + 1; if c_down_right > maxC, c_down_right = 1; end
           potential_neighbors = [potential_neighbors; r+1, c; r+1, c_down_right];
        end
    else
        if r > 1 % Check Up
           c_up_left = c - 1; if c_up_left < 1, c_up_left = maxC; end
           potential_neighbors = [potential_neighbors; r-1, c_up_left; r-1, c];
        end
        if r < maxR % Check Down
           c_down_left = c - 1; if c_down_left < 1, c_down_left = maxC; end
           potential_neighbors = [potential_neighbors; r+1, c_down_left; r+1, c];
        end
    end
    neighbors = potential_neighbors;
end

% ---------------------------------------------------------
% Helper: Connectivity Stats (for Well Mixed)
% ---------------------------------------------------------
function [deg, k_same, k_up, k_down] = get_hex_connectivity(r, total_rows)
    k_same = 2; 
    k_up = 0; if r > 1, k_up = 2; end
    k_down = 0; if r < total_rows, k_down = 2; end
    deg = k_same + k_up + k_down;
end