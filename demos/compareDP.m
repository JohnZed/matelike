function allres = compareDP(nreps, truetheta, NT, usePsi2)
% Performs simple Monte Carlo simulations to compare Empirical Likelihood and
% 2-step GMM estimation of a dynamic panel data model.

% Initialize random seed for repeatable results
randn('state', 1127);

% Use both Psi1 and Psi2 moment conditions
nobs = 1500;
if ~exist('NT')
  NT = 10;
end
if ~exist('usePsi2')
  usePsi2 = true;
end

res_el = [];
res_gmm = [];
ii=1;
while ii <= nreps
  t1 = cputime();
  res = dynamicPanel(nobs, NT, truetheta, 'EL', usePsi2);
  res.time = cputime() - t1;
  if res.flag ~= 0
    disp('Convergence failure');
    continue
  end
  res_el = [res_el ; res ];
  
  t1 = cputime();
  res = dynamicPanel(nobs, NT, truetheta, 'GMM', usePsi2);
  res.time = cputime() - t1;
  res_gmm = [res_gmm ; res ];
  ii = ii+1;
end

fprintf('\nEmpirical likelihood, %d obs\n', nobs);
summarizeres(res_el);

fprintf('\n\nTwo-step GMM, %d obs\n', nobs);
summarizeres(res_gmm);

allres.res_el = res_el;
allres.res_gmm = res_gmm;
allres.nT = NT;
allres.nObs = nobs;

function summarizeres(res)
% Print some stats on the MC experiment
bias = [res.theta]' - repmat(truetheta, nreps, 1);
runtimes = [res.time];
iters = [res.numiter];
estSd = std([res.theta]);

fprintf('Error (absmean/absmedian):   %7.5f  /  %7.5f\n', ...
        mean(abs(bias)), median(abs(bias)));
fprintf('Bias (median):   %7.5f\n', ...
        median(bias));
fprintf('SD of est:    %7.5f\n', estSd);
fprintf('Root mse:     %7.5f\n', sqrt(sum(bias .^ 2) / nreps));
fprintf('Mean/max iterations:  %d  / %d\n', mean(iters), max(iters));
fprintf('Mean time (sec):  %7.5f\n', mean(runtimes));
fprintf('Min/max time (sec):  %7.5f  / %7.5f\n', min(runtimes), max(runtimes));
fprintf('Number of simulations:   %d\n', size(res,1));

end

end