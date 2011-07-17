function res = elSolve(elike, meth, thetaGuess, pguess, varargin)
% Solves the empirical likelihood (or other CR-family) model in elike
%
% Primary usage:
%
%    "res = elSolve(elike, method, thetaGuess, pguess, ...)"
%
% Inputs:
% 
%   "elike" is a structure containing the problem definition. It must
%           be created by:  elike = elSetup(...);  See "help elSetup"
%
%   "method" specifies which objective function to use. It can be 'EL'
%           for empirical likelihood, 'ET' for exponential tilting, or a
%           numerical value indicating the Cressie-Read lambda parmeter to
%           use. If 'method' is 'GMM', a simple two-stage GMM estimator
%           will be used (ignoring pguess and other EL-specific params).
%            
%   "thetaGuess" is the starting value for the theta parameters.
%
%   "pguess" is an optional starting value for the pi parameters. If this
%            argument is empty ([]) or omitted, pi=1/N will be used.
%
%
% Outputs:
%
%   "res" is a structure containing:
%
%      "res.theta"   -- estimated value for theta
%      "res.lagmult" -- lagrange multipliers on moment constraints
%      "res.p"       -- estimated pi for each observation
%      "res.numiter" -- number of iterations taken
%      "res.fval"    -- value of objective function at optimum
%      "res.flag"    -- set to 0 if optimization succceeded, nonzero otherwise
%
% Details:
%
%   elSolve will use the 'zipsolver' matlab package by default. If you have
%   Ipopt and its matlab interface installed, set "elike.solver = 'ipopt'" to
%   use it (Ipopt is more robust and, for large problems, faster than
%   zipsolver).
%
%   See 'help elSetup' for more info on configuring EL problems.
%   
% Author: John Zedlewski (jzedlewski@hbs.edu)  
% 
  global elLoaded;
  if ~(elLoaded)
    error('You must load the matElike package with "elLoad" first');
  end

  nObs = elike.nObs;
  nTheta = elike.nTheta;
  nMom = elike.nMom;
  extraCons = size(elike.linConsLHS, 1);
  nCons = nMom+1 + extraCons;

  % Rescaling to improve stability, auto-computed later
  elObjscale = 1.0;
  elMomscale = 1.0;
  autoscale = false;
  
  % Confirm that our constraints make sense
  if (extraCons && size(elike.linConsLHS,2) ~= nTheta)
    error('Linear constraints must have same size as theta');
  end
  
  % Initial guess
  if nargin < 3 || isempty(thetaGuess)
    theta0 = ones(nTheta, 1);
  else
    printVerbose('Using user-supplied guess for theta\n');
    theta0 = thetaGuess;
  end
  if nargin < 4 || isempty(pguess)
    p0 = ones(nObs, 1) / nObs;
  else
    printVerbose('Using user-supplied guess for p\n');
    p0 = pguess;
  end
  
  % Find the cressie-read lambda corresponding to GEL method
  if isnumeric(meth)
    crLambda = meth;
  else
    % GMM code is in a different function, just call it and return
    if ~isempty(regexpi(meth,'.*GMM'))
      if ~exist('thetaGuess')
        thetaGuess = [];
      end
      switch upper(meth)
        case {'GMM','2GMM'}
          printVerbose('Using 2-step GMM\n');
          elike.gmmType = '2GMM';
        case 'IGMM'
          printVerbose('Using iterated GMM\n');
          elike.gmmType = 'iterated';
        case {'1GMM','ONESTEP'}
          printVerbose('Using 1-step GMM\n');
          elike.gmmType = 'oneStep';
        otherwise
          fprintf('Unknown GMM type %s\n', meth);
          assert(false);
      end
      res = gmmSolve(elike, theta0);
      return;
    end
    
    % Otherwise, we're using a string to represent a GEL method
    switch upper(meth)
      case {'EL','ELIKE'}
        printVerbose('Using empirical likelihood\n');
        crLambda = 0;
      case {'ET','ETILT'}
        printVerbose('Using exponential tilting\n');
        crLambda = -1;
      case 'CUE'
        printVerbose('Using CUE (lambda=-2, neg pi allowed)');        
        crLambda = -2;
      otherwise
        assert(false, 'Unknown method: %s\n', meth);
    end
  end
  
  % Set our objective function based on the lambda parameter
  objFun = @(p,elike,deriv) elFamObj(crLambda, p, elike, deriv);

  % Just evaluate the function and derivatives at the given point
  % if we got the 'eval' option
  if nargin >= 5 && (strcmpi(varargin{1},'eval') || strcmpi(varargin{1},'evalobj'))
    res.f = -elComputeObj(pguess, thetaGuess, elike);
    if strcmpi(varargin{1},'evalobj')
      % only wanted the objective function, so bail out now
      return;
    end
    
    res.fgrad = -elComputeGrad(pguess, thetaGuess, elike);
    res.c = elComputeCons(pguess, thetaGuess, elike);
    res.cJ = elComputeJac(pguess, thetaGuess, 0, elike);
    res.mom = zeros(nObs, nMom);
    for ii=1:nMom
      res.mom(:,ii) = elike.userCompMom(thetaGuess, elike, ii, (ii==1));
    end
    if nargin == 6
      % Got a lambda input, so compute the Hessian of the Lagrangian
      lagmult = varargin{2};
      res.H = -elComputeHess(pguess, thetaGuess, 1.0, lagmult, 0, ...
                             elike);
    else
      res.H = [];
    end
    return
  end
  
  % No bounds on structural parameters, but
  % p constrained between 0 and 1
  lb = { zeros(nObs,1) ; repmat(-Inf, nTheta, 1) };
  ub = { ones(nObs,1) ; repmat(+Inf, nTheta, 1) };
  fixvars = [];
  if nCons == nTheta
    % Exactly-identified case: p must be 1 / N
    lb(1:nObs) = 1 / nObs;
    ub(1:nObs) = 1 / nObs;
  end
  
  % If we got the 'profile' option, compute profile empirical likelihood 
  % for given theta (holding theta fixed, maximize over p, then evaluate)
  % Requires an extra argument: a logical vector indicating which
  % components of theta to hold fixed for the profile
  if nargin >= 5 && strcmpi(varargin{1},'profile')
    assert(nargin > 5, 'Need to specify theta vars to hold fixed for profile');
    assert(~isempty(thetaGuess), 'Need to specify theta value for profile');
    
    pvin = varargin{2} == 1;
    assert(numel(pvin) == nTheta, 'profvars should be theta-sized vector');
    if strcmpi(elike.solver,'ipopt')
      proflb = repmat(-Inf,nTheta,1);
      profub = repmat(+Inf,nTheta,1);
      proflb(pvin) = thetaGuess(pvin);
      profub(pvin) = thetaGuess(pvin);
      lb{2} = proflb;
      ub{2} = profub;
    else
      assert(false, 'Only Ipopt supports fixed variables currently\n');
    end
  end

  % Constraints set equal to 0
  lbc = zeros(nCons, 1);
  ubc = zeros(nCons, 1);
  
  % Debugging stuff
  prevTheta = [];
  prevF = NaN;
  iterFunc = '';
  if (elike.printTheta)
    iterFunc = @elIterFunc;
  end

  % Rescale the objective function or moments to make
  % the gradient and constrants of order 1
  if ~autoscale
    [grad0 gradt] = elComputeGrad(p0, theta0, elike);
    maxg0 = max(abs(grad0));
    s = 10.0 / maxg0;
    elObjscale = 10 ^ round( log10(s) ); % Round to nearest power of 10
  
    c0 = elComputeCons(p0, theta0, elike);
    s = 10.0 / max(abs(c0));
    elMomscale = 10 ^ round( log10(s) ); % Round to nearest power of 10
  end
  
  % Actually do the optimization with either Ipopt or Zipsolver
  if strcmpi(elike.solver,'ipopt')
    printLev = 1;
    if (elike.verbose)
      printLev = 5;
    end
    if autoscale
      ipscale = 'gradient-based';
    else
      ipscale = 'none';
    end
    [pstar thetastar lagmult numiter] = ipopt({p0,theta0},lb,ub,lbc,ubc,...
                                              @elComputeObj, @elComputeGrad, ...
                                              @elComputeCons, @elComputeJac, @elComputeHess, ...
                                              elike, '', [], ...
                                              'print_level', printLev, ...
                                              'linear_solver', elike.linearSolver, ...
                                              'nlp_scaling_method', ipscale, ...
                                              'mu_strategy', 'adaptive', ...
                                              'max_iter', elike.maxIter);
    lagmult = lagmult.lambda';
    
    % Assume success unless maximum iterations reached  
    res.flag = (numiter == elike.maxIter);
  elseif strcmpi(elike.solver,'fmincon')
    % Use the fmincon function from Matlab's toolbox (needs R2008a)
    if verLessThan('optim','4.0')
      printVerbose(['Optimization toolbox 4.0 (Matlab R2008a) required for ' ...
                    'fmincon solver\n']);
    end
    
    if elike.verbose
      displayopt = 'iter';
    else
      displayopt = 'notify';
    end

    % Must use the interior-point algorithm and provide Hessian
    fminOpts = optimset('Algorithm','interior-point',...
                        'GradObj','on',...
                        'Diagnostics','off',...
                        'Display',displayopt,...
                        'Hessian','user-supplied',...
                        'HessFcn',@elHessFmincon,...
                        'GradConstr','on');
    fminOpts = addExtraOpts(fminOpts,elike.solverOpts);
    
    lb = cell2mat(lb); ub = cell2mat(ub);
    x0 = [p0; theta0];    
    [xstar,fval,flag,output,lambdaout] = fmincon(@elObjZipsolver,x0,...
                                                 [],[],[],[],lb,ub,...
                                                 @elConsFmincon, fminOpts);
    lagmult = lambdaout.eqnonlin;
    pstar = xstar(1:nObs);
    thetastar = xstar(nObs+1:end);
    numiter = output.iterations;
    res.flag = flag;
  elseif strcmpi(elike.solver,'zipsolver')
    opt = zipopts('autoscale', autoscale, 'maxIter', elike.maxIter, ...
                  'printstep', elike.verbose, 'maxWatchfail', 10);
    opt = addExtraOpts(opt, elike.solverOpts);

    if isempty(elike.lambdaGuess)
      lamGuess = 0.01 * ones(nCons,1);
    else
      lamGuess = elike.lambdaGuess;
    end
      
    if (elike.pPositive)
      posvars = [ ones(nObs,1) ; zeros(nTheta,1) ];
    else
      posvars = zeros(nObs + nTheta,1);
    end
    
    x0 = [p0; theta0];
    [xstar,info] = zipsolver(@elObjZipsolver, @elConsZipsolver, @elHessZipsolver, ...
                             x0, opt, lamGuess, posvars);
    pstar = xstar(1:nObs);
    thetastar = xstar(nObs+1:end);
    lagmult = info.lagmult;
    numiter = info.numiter;
    res.flag = info.flag;
  else
    assert(false, 'Unknown solver %s, should be zipsolver, ipopt, or fmincon', ...
           elike.solver)
  end

  res.p = pstar;
  res.theta = thetastar;
  res.pGuess = p0;
  res.thetaGuess = theta0;
  res.lagmult = lagmult;
  res.numiter = numiter;
  res.fval = -elComputeObj(res.p, res.theta, elike) / elObjscale;
  res.elike = elike;
  res.meth = meth;

function f = elComputeObj(p, theta, elike)
%
% Compute the value of the objective function
%
  f = elObjscale * objFun(p, elike, 0);
  
  prevTheta = theta;
  prevF = f;

end;

function [gradp, gradtheta] = elComputeGrad(p, theta, elike)
%
% Compute the gradient of the objective function
%
  gradp = elObjscale * objFun(p, elike, 1);
  gradtheta = elObjscale * zeros(length(theta),1);
end;

function cons = elComputeCons(p, theta, elike)
%
% Compute the constraints based on the moment conditions
%
  cons = zeros(nCons,1);

  % Problem-specific moment conditions
  for ii=1:elike.nMom
    Mii = elike.userCompMom(theta, elike, ii, (ii==1));
    % Check that the matrix looks right
    if (size(Mii,1) ~= elike.nObs || size(Mii,2) ~= 1)
      errmsg = sprintf('Expected moment matrix to have dimension %d x 1\n', ...
                       elike.nObs);
      error(errmsg);
    end
    cons(ii) = p' * Mii;
  end

  % Make p's sum to one
  cons(nMom+1) = 1 - sum(p);
  
  % Any additional linear constrains
  if (extraCons)
    cons(nMom+2:end) = elike.linConsLHS * theta - elike.linConsRHS;
  end
  
  cons = cons * elMomscale;
end;

function H = elComputeHess(p, theta, objscale, lambda, structOnly, elike)
%
% Compute the Hessian of the Lagrangian, using automatic
% differentiation to find the parts that depend on the moment
% conditions, but using the simple structure of the EL objective
% function for the rest of the sparse matrix.
% 
  NP = elike.nParams;
  nMom = elike.nMom;
  nObs = elike.nObs;
  objscale = objscale * elObjscale;

  if structOnly
    Hupper = spdiags( ones(nObs,1), 0, nObs, NP );
    Hlower = [ ones(nTheta, nObs), tril(ones(nTheta,nTheta)) ];
    H = sparse( [ Hupper ; Hlower ] );
  else
    % Second deriv of obj func wrt p, scaled by objscale
    pdiag = objscale * objFun(p, elike, 2);
    
    Hupper = spdiags( pdiag , 0, nObs, NP );
    
    % Compute constraint derivs and cross-partials
    hthetaP = zeros(nTheta, nObs);     % theta-p cross-partial
    htheta2d = zeros(nTheta,nTheta);       % theta-theta second deriv
    
    theta_ad = hessianinit(theta);
    for ii=1:nMom
      if (isempty(elike.userCompHessTheta) || isempty(elike.userCompJac))
        % AD approach
        M_ii = elike.userCompMom(theta_ad, elike, ii, (ii==1));
        Mp = p' * M_ii;
        
        j_ii = M_ii.dx;
        hess_ii = Mp.hx;
      else
        % Custom function approach
        M_ii = elike.userCompMom(theta, elike, ii, (ii==1));
        Mp = p' * M_ii;
        j_ii = elike.userCompJac(theta, elike, ii, (ii==1));
        hess_ii = elike.userCompHessTheta(theta, elike, ii, (ii==1));
      end
      
      htheta2d = htheta2d + lambda(ii)  * hess_ii;
      hthetaP = hthetaP + lambda(ii) * (j_ii');
    end
    
    Hlower = elMomscale * sparse([hthetaP tril(htheta2d)]);
    H = [ Hupper ; Hlower ];
  end
end;

function J = elComputeJac(p, theta, structOnly, elike)
%
% Computes the Jacobian of the constraints
%
  nObs = elike.nObs;
  nMom = elike.nMom;

  if structOnly
    if (isempty(elike.momSparsePat))
      % Assume dense jacobian by default
      J = ones(nCons, elike.nParams);
    else
      % Allow user to override with own sparsity pattern
      J = [ elike.momSparsePat', ones(nMom,nTheta) ; 
            ones(1,nObs), zeros(1, nTheta) ];
    end
  else
    if isempty(elike.momSparsePat)
      J = zeros(nCons, elike.nParams);
    else
      J = [ elike.momSparsePat', ones(nMom,nTheta) ; 
            ones(1,nObs), zeros(1, nTheta) ];
    end
    
    % Fill in the blocks of the Jacobian for each moment
    theta_ad = gradientinit(theta);
    for ii=1:nMom
      if isempty(elike.userCompJac)
        M_ad_ii = elike.userCompMom(theta_ad, elike, ii, (ii==1));
        j_ii = M_ad_ii.dx;
        m_ii = M_ad_ii.x;
      else
        j_ii = elike.userCompJac(theta, elike, ii, (ii==1));
        m_ii = elike.userCompMom(theta, elike, ii, (ii==1));
      end
      
      % Deriv wrt theta
      if nTheta > 0
        J(ii,nObs+1:end) = j_ii' * p;
      end
      
      % Deriv wrt p
      J(ii,1:nObs) = m_ii';
    end

    % Deriv of condition that p's add up to 1
    J(nMom+1,1:nObs) = -1;
    
    % Additional linear constraints on theta(if any)
    if (extraCons)
      J(nMom+2:end,nObs+1:end) = elike.linConsLHS;
    end
  end

  J = elMomscale * sparse(J);
end;


function elIterFunc(repnum, fval, elike)
% Print the most recently used theta values
  if (prevF ~= fval)
    fprintf('prevF=%g, fval=%g, printed theta may be incorrect\n', ...
            prevF, fval);
  end
  fprintf('theta[%3d]:  ', repnum);
  fprintf('%9.6f  ', prevTheta);
  fprintf('\n');
end


%%
%% Functions to bridge between the ipopt interface and the zipsolver
%% interface (they just call the relevant elComputeXXX functions) 
%%

function [f,g] = elObjZipsolver(x)
  p = x(1:nObs);
  theta = x(nObs+1:end);
  f = elComputeObj(p, theta, elike);
  if nargout > 1
    [gp,gt] = elComputeGrad(p, theta, elike);
    g = [gp ; gt];
  end
end; % end function elObjZipsolver

function H = elHessZipsolver(x, objscale, lambda)
  p = x(1:nObs);
  theta = x(nObs+1:end);
  H = elComputeHess(p, theta, objscale, lambda, false, elike);
  
  % Convert from a lower triangle to a full matrix
  Htri1 = tril(H,-1);
  Htri = tril(H);
  H = Htri + Htri1';
end; % end function elHessZipsolver

function [c,A] = elConsZipsolver(x)
  p = x(1:nObs);
  theta = x(nObs+1:end);
  c = elComputeCons(p, theta, elike);
  if nargout > 1
    A = elComputeJac(p, theta, false, elike);
  end
end; % end function elConsZipsolver

function opt = addExtraOpts(opt, newOpts)
% Merge name-value pair options from 'newOpts' into opt
  soptfields = fields(newOpts);
  for iSol=1:length(soptfields)
    solopt = soptfields{iSol};
    optval = newOpts.(solopt);
    opt.(solopt) = optval;
    printVerbose('Set solver option %s to : ', solopt);
    if (elike.verbose)
      disp(optval)
    end
  end  
end

function printVerbose(fmt, varargin)
  if elike.verbose
    fprintf(fmt, varargin{:});
  end
end


%% 
%% Bridge functions for fmincon (very similar to zipsolver)
%%
function [cIneq,cEq,cJIneq,cJEq] = elConsFmincon(x)
  cIneq = [];
  cJIneq = [];
  
  [cEq,cJEq] = elConsZipsolver(x);
  cJEq = cJEq';
end

function H = elHessFmincon(x, lambdaIn)
  lambda = lambdaIn.eqnonlin;
  H = elHessZipsolver(x, 1.0, lambda);
end


end