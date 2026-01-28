Nreps = 10000;
results = zeros(Nreps,2);
numCells = 16;
lambda = 0.5;

tic;
parfor trial = 1:Nreps
     [w, t] = MutantSimulation(numCells, 1, lambda);
     results(trial, :) = [w, t]
end
toc;
probToMTDominate = sum(results(:,1) == -1) / Nreps;
disp(probToMTDominate)

% Mutant wins
figure;
subplot(1,2,1)
mask = results(:,1) == -1;
histogram(results(mask, 2), 50);
xlabel('Time to Mutant Domination');
ylabel('Frequency');
title(sprintf('Distribution of Monoclonal Times for %d cells (N=10000)', numCells));
subplot(1,2,2)
sorted_steps = sort(results(mask,2));
cumulative_percent = (1:numel(sorted_steps)) ./ numel(sorted_steps) * 100;
plot(sorted_steps, cumulative_percent, '-k', 'LineWidth', 1.5);
ylim([0,100])
xlabel('Time')
ylabel('Monoclonal Crypt (%)')

% Wild Type wins
figure;
subplot(1,2,1)
mask = results(:,1) ~= -1;
histogram(results(mask, 2), 50);
xlabel('Time to Wild Type Domination');
ylabel('Frequency');
title(sprintf('Distribution of Monoclonal Times for %d cells (N=10000)', numCells));
subplot(1,2,2)
sorted_steps = sort(results(mask,2));
cumulative_percent = (1:numel(sorted_steps)) ./ numel(sorted_steps) * 100;
plot(sorted_steps, cumulative_percent, '-k', 'LineWidth', 1.5);
ylim([0,100])
xlabel('Time')
ylabel('Monoclonal Crypt (%)')

figure;
histogram(results(:, 1), -1:numCells + 1);
xlabel('Winning Cell Label (Original Position)');
ylabel('Frequency');
title('Win Count per Initial Position');


