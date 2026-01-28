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
figure;
% Simulated points with error bars (SEM)
errorbar(lambdas, pDom, errLow, errHigh, 'o', 'LineWidth', 1.2, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'auto', 'DisplayName', 'Simulation (95% CI)');
% plot(lambdas, pDom, '-o', 'LineWidth', 1.5);
hold on;
fplot(theorectic_function, [0,3.5])
xlabel('Lambda');
ylabel('Prob(mutant dominates)');
title(sprintf('Mutant Win Probability (N = %d, N_{reps} = %d)', numCells, Nreps));
legend('Simulated probability', 'Theoretic function of probability')
xlim([-0.3,4])
ylim([-0.1,0.8])
grid on;
