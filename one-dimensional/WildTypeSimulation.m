function [timeMonoclonal, winningCell] = WildTypeSimulation(visualize, numCells)
% Run the simulation on a base vector representing the base of the crypt.
% With the strategy of picking one cell uniformly to spread and replace
% either the neighbouring left or right cell.

% -- Parameters --
N = numCells;
cells = 1:N;
maxIter = 1e6;

history = cells;
time_step = 0;

% Continue until Monoclonal Conversion
while length(unique(cells)) > 1 && time_step < maxIter
    time_step = time_step + 1;
    
    % Pick a cell to divide
    dividing_idx = randi(N);
    label_to_spread = cells(dividing_idx);
    
    % 0 for replace left neighbor, 1 for replace right neighbor
    direction = randi([0, 1]); 
    
    if direction == 0 
        % Replace Left Neighbor
        neighbor_idx = dividing_idx - 1;
        if neighbor_idx < 1
            neighbor_idx = N; % Periodic boundary
        end
    else 
        % Replace Right Neighbor
        neighbor_idx = dividing_idx + 1;
        if neighbor_idx > N
            neighbor_idx = 1; % Periodic boundary
        end
    end

    cells(neighbor_idx) = label_to_spread;
    history = [history; cells];
end
if time_step >= maxIter
    error('Exceeded max iterations');
end
winningCell = cells(1);
timeMonoclonal = time_step;

if visualize == true
    disp(history)
    fprintf('Monoclonal conversion reached at step: %d\n', time_step);
    fprintf('Winning clone: %d\n', cells(1));
end

end

