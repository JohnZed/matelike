function res = gmmSolve(elike, thetaGuess, varargin)
% Solves for theta using 2-step GMM.
%
% Usage:
%
%   res = gmmSolve(elike, thetaGuess, ...)
%
% Inputs:
%
%   'elike' is a structure containing the problem definition. It must
%         be created by:  elike = elSetup(...);  See 'help elSetup'
%         The initial weighting matrix may be set by assigning
%         it to 'elike.W1'. If it is not set, the identity matrix is used.
%
%   'thetaGuess' is an optional starting guess for theta.
%
% Outputs:
%
%   'res' is a structure like that returned by elSolve(...), containing
%         'res.theta' with the estimated value of theta.
%
% Details:
%
%    This is not a very sophisticated GMM implementation and is included only
%    for simple comparisons. The second-step weighting matrix is stored in
%    'res.elike.W'.
%

nObs = elike.nObs;
nTheta = elike.nTheta;
nMom = elike.nMom;

% No bounds on structural parameters by default
lb = { repmat(-Inf, nTheta, 1) };
ub = { repmat(+Inf, nTheta, 1) };

% No constraints (GMM)
lbc = [];
ubc = [];

% Initial guess
if nargin < 2 || isempty(thetaGuess)
  theta0 = ones(nTheta, 1);
else
  printVerbose('Using user-supplied guess for theta\n');
  theta0 = thetaGuess;
end

extraargs = varargin;
if elike.zeroHess
  extraargs = { extraargs{:}, 'hessian_constant', 'yes' };
end

% First step weighting matrix
W1 = fieldopt(elike, 'W1', eye(nMom));
printVerbose('GMM first stage');

% Stage 1 optimization
[thetastar1 lagmult1 numiter1] = doGmmSolve(theta0, W1);

res.flag = 0;

if fieldoptstr(elike, 'gmmType', 'oneStep')
  % Only do a single step
  printVerbose('Only one step GMM\n');
  W2 = W1;
  thetastar = thetastar1;
  lagmult = lagmult1;
  numiter = numiter1;
  res.gmmIter = 1;
elseif fieldoptstr(elike, 'gmmType', 'iterated')
  % Iterate the weighting matrix until convergence
  maxGmmIter = fieldopt(elike, 'maxGmmIter', 100);
  convergeTol = fieldopt(elike, 'covergeTol', 1e-4);
  prevTheta = thetastar1;
  for gmmIter=1:maxGmmIter
    W = gmmComputeWeights(prevTheta);
    [thetastar lagmult numiter] = doGmmSolve(prevTheta, W);
    tchange = max(abs(thetastar - prevTheta));
    printVerbose('Iteration %3d theta change: %9.6f\n', gmmIter+1, tchange);
    if tchange < convergeTol
      % Convergence achieved
      break
    end
    prevTheta = thetastar;
  end
  if gmmIter == maxGmmIter
    fprintf('Convergence not achieved after %d GMM iterations, recent change %9.5f\n', ...
      maxGmmIter, tchange);
    res.flag = -1;
  end
  res.gmmIter = gmmIter;
  W2 = W;
else
  % By default, use two-step GMM
  W2 = gmmComputeWeights(thetastar1);
  printVerbose('Second stage optimization\n');
  [thetastar lagmult numiter] = doGmmSolve(thetastar1, W2);
  res.gmmIter = 2;
end

res.theta = thetastar;
res.lagmult = lagmult;
res.numiter = numiter;
res.thetaGuess = theta0;
res.elike = elike;
res.meth = 'GMM';
Wstar = gmmComputeWeights(thetastar);
res.V = gmmComputeVar(thetastar, Wstar, elike);
res.W1 = W1;
res.W2 = W2;
res.fval = gmmComputeObj(thetastar, elike);

  function [t,l,n] = doGmmSolve(ts, W)
    elike.W = W;
    if (strcmpi(elike.solver,'ipopt'))
      printLev = 1;
      if (elike.verbose)
        printLev = 5;
      end
      [t l n] = ipopt(ts,lb,ub,lbc,ubc,...
                      @gmmComputeObj, @gmmComputeGrad, ...
                      '', '', @gmmComputeHess, ...
                      elike, '', [], ...
                      'mu_strategy', 'adaptive', ...
                      'max_iter', elike.maxIter, ...
                      'print_level', printLev, ...
                      'linear_solver', elike.linearSolver, ...
                      extraargs{:});
    else
      opt = zipopts(varargin{:});
      opt.printstep = elike.verbose;
      posvars = [];
      [xstar, info] = zipsolver(@gmmObjZipsolver, [], ...
        @gmmHessZipsolver, ts, opt, ...
        [],[], elike);
      t = xstar;
      l = info.lagmult;
      n = info.numiter;
    end
  end

  function f = gmmComputeObj(theta, elike)
    m = [];
    for ii=1:nMom
      m_obs = elike.userCompMom(theta, elike, ii, (ii==1));
      mii = sum(m_obs) / elike.nObs;
      m = [m ; mii];
    end
    
    f = m' * elike.W * m;
  end


  function gradtheta = gmmComputeGrad(theta, elike)
    theta_ad = gradientinit(theta);
    f_ad = gmmComputeObj(theta_ad, elike);
    gradtheta = f_ad.dx;
  end

  function H = gmmComputeHess(theta, objscale, lambda, structOnly, ...
      elike)
    if structOnly
      % Use a full matrix
      H = sparse( tril(ones(nTheta,nTheta)) );
    else
      theta_ad = hessianinit(theta);
      f_ad = gmmComputeObj(theta_ad, elike);
      H = sparse( tril(f_ad.hx) );
    end

  end

  function W = gmmComputeWeights(theta)
    % Update the weighting matrix with the gradients
    m1 = zeros(nObs, nMom);
    for ii=1:nMom
      m1(:,ii) = elike.userCompMom(theta, elike, ii, (ii==1));
    end
    opgrad = m1' * m1 / elike.nObs;
    W = inv(opgrad);
  end

  function V = gmmComputeVar(theta, W, elike)
    theta_ad = gradientinit(theta);
    m = [];
    for ii=1:nMom
      m_obs = elike.userCompMom(theta_ad, elike, ii, (ii==1));
      m = [m , m_obs];
    end
    mbar = mean(m);
    G = mbar(:).dx;
    mx = m.x;
    S = mx' * mx / nObs;

    bread = inv(G' * W * G);
    meat = (G' * W) * S * (W * G);
    V = bread * meat * bread / nObs;
  end

%%%%% ZIPsolver support %%%%%

  function H = gmmHessZipsolver(theta, objscale, lambda, elike)
    H = gmmComputeHess(theta, objscale, lambda, false, elike);

    % Convert from a lower triangle to a full matrix
    Htri1 = tril(H,-1);
    Htri = tril(H);
    H = Htri + Htri1';
  end


  function [f,g] = gmmObjZipsolver(theta, elike)
    f = gmmComputeObj(theta, elike);
    if nargout > 1
      g = gmmComputeGrad(theta, elike)';
    end
  end


% Utility functions

  function printVerbose(fmt, varargin)
    if elike.verbose
      fprintf(fmt, varargin{:});
    end
  end

end

