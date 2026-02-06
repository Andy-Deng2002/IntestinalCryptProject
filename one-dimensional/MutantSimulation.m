function [winningCell, timeToMonoclonal] = MutantSimulation(numCells, MutantPos, lambda, visualize)

if nargin < 4
    visualize = false; % Set default visualization to false
end
if nargin < 3
    MutantPos = randi(numCells); % Default: random select one cell to be mutant
end

if MutantPos > numCells || MutantPos < 1
    error("Mutant Position should be within the cells array")
end

% -- Parameters --
N = numCells;
cells = 1:N;
maxIter = 1e6;

% Initialize the cells array with mutant cell included
cells(MutantPos) = -1;

history = cells;
time_step = 0;

% Continue until Monoclonal Conversion
while length(unique(cells)) > 1 && time_step < maxIter
    time_step = time_step + 1;

    % Set up the division distribution
    weights = (cells == -1) * (lambda) + (cells ~= -1) * 1;
    weights = max(weights, 0);
    dist = weights / sum(weights); 
    
    % Pick a cell to proliferate
    proliferate_idx = randsample(N, 1, true, dist);
    
    if rand < 0.5 
        % occupy left neighbour
        removed_idx = proliferate_idx - 1;
        if removed_idx < 1, removed_idx = N; end
    else 
        % occupy right neighbour
        removed_idx = proliferate_idx + 1;
        if removed_idx > N, removed_idx = 1; end
    end

    cells(removed_idx) = cells(proliferate_idx);
    if visualize
        history = [history; cells];
    end
end
if time_step >= maxIter
    error('Exceeded max iterations');
else
    winningCell = cells(1);
end
timeToMonoclonal = time_step;

if visualize == true
    disp(history)
    fprintf('Monoclonal conversion reached at step: %d\n', time_step);
    fprintf('Winning clone: %d\n', cells(1));
end

end