% demoPoissinst --
%
%   Estimates an IV-Poisson model of cigarette consumption
%

% Load the Mullahy (1997) dataset, which contains X, y, and Z
% This dataset is derived from 'ivcig', provided with the Stata module
% 'ivpois' by Austin Nichols (2007).

load ivcig;

nObs = size(X,1)
keepers = ((1:nObs) < 1000)';  % shrink this number to use a subset
X = X(keepers,:);
y = y(keepers,:);
Z = Z(keepers,:);
nObs = size(X,1)

% Add a constant term
X = [X, ones(nObs,1)];

% Optional: Interact income and education with various things
% zInc = repmat(X(:,4),1,6) .* X(:,[2,3,5,7,9,10]) / 1000;
% zEduc = repmat(X(:,7),1,4) .* X(:,[2,3,9,10]) / 100;
zInc = [];
zEduc = [];

% Only the first element of X (k210) is endogenous
Z = [Z, zInc, zEduc, X(:,2:end)];

% Actually solve the model two different ways (see poissinst.m for details)
fprintf('Num X: %d     Num inst:  %d\n', size(X,2), size(Z,2));
tic
  resGmm = poissinst(X, y, Z, 'GMM');
toc
tic
  resEl = poissinst(X, y, Z, 'EL');
toc

% View the results from the two methods
names = strvcat(xnames{:}, 'intercept');
elStderr = elModelSumm(resEl, names);
gmmStderr = elModelSumm(resGmm, names);

% How much do they differ?
relDiff = abs(resGmm.theta - resEl.theta) ./ elStderr

% Print as a LaTeX table
fullnames = {'Habit stock', 'Price 79', 'Rest res. 79', 'Income', 'Age', ...
             '$Age^2$/1000', 'Educ', '$Educ^2$', 'Fam size', 'White', 'Intercept'};
mres = [ resEl.theta, elStderr, resGmm.theta, gmmStderr, relDiff];
for ii=1:size(mres,1)
  fprintf('%s & %9.4f & (%.4f) & %9.4f & (%.4f) & %9.4f \\\\ \n', ...
          fullnames{ii}, mres(ii,:));
end
