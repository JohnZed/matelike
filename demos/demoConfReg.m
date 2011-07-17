% Demonstrates computing and plotting confidence intervals.

% Solve a linear model and store the result as 'res'
runSimpleDemo
res.elike.verbose = false;

% Default: compute 95% confidence region by grid search
tic
  [cintL,cintU] = elConfRegion(res);
toc
fprintf('95%% confidence interval: \n');
for ii=1:length(cintL)
  fprintf('    [%8.4f, %8.4f]\n', [cintL(ii),cintU(ii)]);
end

% Plot the profile of p-values with both 90% and 99% regions
opts.varnames = { '\theta_1', '\theta_2' };
opts.rotateLabels = true;  % Fit the text better on the x axis
opts.plotPoints = true;    % See which points were sampled exactly
opts.labelCI = true;
[cintL2,cintU2] = elConfRegion(res, [1,2], [0.90,0.99], true, opts);


