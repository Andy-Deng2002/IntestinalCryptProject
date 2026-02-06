numCells = 16;
MutantPos = 8;
Nreps = 4000;
lambdas = linspace(0.001,3.5,25);
pDom = zeros(size(lambdas));
ciLow = zeros(size(lambdas));
ciHigh = zeros(size(lambdas));

tic;
z = 1.96;
parfor i = 1:numel(lambdas)
    lam = lambdas(i);
    wins = false(Nreps,1);
    results = zeros(Nreps,2);
    for trial = 1:Nreps
        [w, t] = MutantSimulation(numCells, MutantPos, lam, false);
        results(trial,:) = [w,t];
        wins(trial) = (w == -1);
    end

    valid = ~isnan(results(:,1));
    n = nnz(valid);
    k = sum(wins(valid));
    p = k / n;
    se = sqrt(p*(1-p)/n);
    lo = p - z * se;
    hi = p + z * se;
    ciLow(i) = lo;
    ciHigh(i) = hi;
    pDom(i) = p;
end
toc;

errLow = pDom - ciLow;
errHigh = ciHigh - pDom;


theorectic_function = @(lambda) (1 - lambda^(-1)) / (1-lambda^(-numCells));

% Create a professional looking figure
figure('Color', 'w', 'Position', [100, 100, 800, 550]);
hold on;

% Plot theoretical function (smooth line)
% Using a distinct orange color
hLine = fplot(theorectic_function, [0, 3.5], 'LineWidth', 2, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', 'Theoretic function of fixation probability');

% Plot simulated points with error bars
% Using filled blue markers
hData = errorbar(lambdas, pDom, errLow, errHigh, 'o', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerEdgeColor', [0 0.4470 0.7410], ...
    'Color', [0 0.4470 0.7410], ... % Line color for error bars
    'CapSize', 8, ...
    'DisplayName', 'Simulated probability (95% CI)');

xlabel('Proliferative Advantage (\lambda)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability of Mutant Fixation', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Mutant Fixation Probability\n(N = %d, N_{reps} = %d)', numCells, Nreps), ...
    'FontSize', 14);

legend([hData, hLine], 'Location', 'southeast', 'FontSize', 11);

% Tighter limits and nicer grid
xlim([0, 3.6]);
ylim([-0.05, 0.85]); 
grid on;
set(gca, 'Box', 'on', 'FontSize', 12, 'LineWidth', 1.1);
hold off;
saveFigure