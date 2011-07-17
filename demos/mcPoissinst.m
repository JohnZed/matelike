function allres = mcPoissinst(nReps, method)
% Perform a monte carlo simulation to evaluate performance
% of Poisson-IV estimator
  nObs = 1000;
  sigX = 1;
  nXendog = 1;
  nXexog = 10;
  nX = nXexog + nXendog;
  nZ = 10;
  truebeta = 0.5 * (1:nX)' - 1;

  randn('state', 9999);
  rand('state', 9999);

  allres = zeros(nReps,nX);
  for iRep=1:nReps
    % We have exogenous instruments Z
    Zextra = [ sigX * randn(nObs,nZ - 1), ones(nObs,1)];
    
    % Xendog is a function of Z plus error
    sigXeps = 0.75;
    thetaXZ = (0.10 / nZ) * (1:nZ)' * (1:nXendog);
    xeps = sigXeps * randn(nObs, nXendog);
    Xendog = Zextra * thetaXZ + xeps;
    
    % Form the full structural and instrument matrices
    Xexog = [sigX * randn(nObs, nXexog-1) ];
    X = [Xendog, Xexog, ones(nObs, 1) ];
    Z = [ Zextra, Xexog ];
    
    % Eta is an unobservable correlated with Xendog
    logeta = xeps * 0.5 * ones(nXendog,1);
    
    % Dependent variable is a Poisson count, based on
    % an exponential mean in X and the unobservable eta
    y = poissrnd( exp(X * truebeta + logeta) );

    % Actually solve the model
    res = poissinst(X, y, Z, method, false);
    
    allres(iRep,:) = res.theta;
  end

  resDiff = allres - repmat(truebeta',nReps,1);
  fprintf('Method:  %s\n', method);
  fprintf('True beta:    '); disp(truebeta')
  fprintf('Mean diff:    '); disp(mean(resDiff))
  fprintf('Median bias:  '); disp(median(resDiff))
  fprintf('Root MSE:     '); disp(sqrt(sum(resDiff .^ 2) / nReps))
  fprintf('Std dev:      '); disp(std(allres))

end