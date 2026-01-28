Nreps = 10000;
results = zeros(Nreps,2);
numCells = 16;

tic;
parfor trial = 1:Nreps
    [time, WinningCell] = WildTypeSimulation(false, numCells);
    results(trial, :) = [time, WinningCell]; 
end
toc;

figure;
histogram(results(:, 1), 50);
xlabel('Time Steps to Monoclonal');
ylabel('Frequency');
title(sprintf('Distribution of Monoclonal Times for %d cells (N=10000)', numCells));

figure;
histogram(results(:, 2), 1:numCells + 1);
xlabel('Winning Cell Label (Original Position)');
ylabel('Frequency');
title('Win Count per Initial Position');

figure;
sorted_steps = sort(results(:,1));
cumulative_percent = (1:Nreps) ./ Nreps * 100;
plot(sorted_steps, cumulative_percent, '-k', 'LineWidth', 1.5);
ylim([0,100])
xlabel('Time')
ylabel('Monoclonal Crypt (%)')
