function elike = elSetup(nObs, nMom, nTheta, userCompMom, varargin)
% Initializes an EL (or other Cressie-Read) problem
%
% Usage:
%
%   "elike = elSetup(nObs, nMom, nTheta, momFunction, ...)"
%
% Inputs:
%     "nObs, nMom, nTheta" are the number of observations, number of
%         moment conditions, and number of elements in theta
%
%     Extra arguments may follow 'momFunction' to set various options
%         for the optimization process. Options are specified as pairs,
%         where the first argument is the option name and the second
%         is the option value, like:
%
%             "elSetup(nobs,nmom,ntheta,func, 'verbose', false)"
%
%         which sets the 'verbose' option to false, greatly reducing the
%         output displayed.
%
%     "momFunction" is a function handle to compute the moment
%         conditions for each observation, with the form:
%      
%        "M = momFunction(theta, elike, mnum, newTheta)"
%
%          "mnum" indicates which moment to evaluate, "theta" is the point
%          at which to evaluate it, and "elike" is the same structure
%          that was returned from elSetup. If you need to pass additional
%          paraneters to your function, you can stash them in "elike" as
%          fields.
%
%          "newTheta" is true if this theta differs from the one
%          that was passed in the most recent call. If you need to
%          precompute some values only once for each theta, "newTheta"
%          tells you when it's time to redo that computation.
%
%          The returned value M should be an nObs x 1 vector containing
%          moment condition "mnum" evaluated at each observation.
%          
%     *** IMPORTANT restrictions on momFunction ***
%
%          elSolve will use automatic differentiation to compute the
%          derivatives of momFunction. For smooth functions computed
%          analytically, this should work transparently. It will NOT work if
%          momFunction uses some iterative algorithm internally, is
%          non-differentiable, or uses some obscure external functions that
%          are not understood by the differentiation library (Intlab).
%
%          The matrix returned by momFunction must be of the special AD
%          type. This will happen by automatically if the result is produced by
%          operations involving "theta", but it will fail if you manually
%          initialize M with something like "M = zeros(nObs,1)".  If you find
%          this initialization convenient, you can use:
%             "M = zeros(nObs,1) * theta(1)",
%          which will ensure that M has the correct type.
%
%          If automatic differentiation is not suitable for your model, you
%          can provide functions in "elike.userCompJac" and
%          "elike.userCompHessTheta" to compute them. (See elLinearIV for
%          an example.)
%
% Output:
%
%     "elike" is a structure that can be passed to elSolve
%
% Requirements:
%
%     You must have Intlab installed to use elSetup. You also need
%     to call "elLoad" sometime before elSetup to initialize the
%     matelike library.
%
% Author:  John Zedlewski (jzedlewski@hbs.edu)
%

if verLessThan('matlab','7.4')
  fprintf('Warning: matElike not tested with your version of Matlab\n');
end

elike.nObs = nObs;
elike.nMom = nMom;
elike.nTheta = nTheta;
elike.nParams = nObs + nTheta;

% Function handles
elike.userCompMom = userCompMom;
elike.userCompJac = [];
elike.userCompHessTheta = [];

% Linear constraints: linConsLHS * theta = linConsRHS
elike.linConsLHS = [];
elike.linConsRHS = [];

% Special: constrained objective function, minimize linear function of theta
elike.objCons = [];
elike.thetaMinFunc = [];

% Set to false if you do not want to constrain p to be >= 0
elike.pPositive = true;
    
% Set verbose to false to hide output
elike.verbose = true;

% Misc configuration
elike.printTheta = false;
elike.zeroHess = false;
elike.maxIter = 50;
elike.lambdaGuess = [];

% Fill this in if the moments are sparse with a known pattern
elike.momSparsePat = [];

% Must have ipopt mex interface installed if you
% change solver to 'ipopt'
global elDefaultSolver;
elike.solver = elDefaultSolver;

% linearSolver only applies to Ipopt. May set to 'pardiso' if that solver
% is installed.
elike.linearSolver = 'ma27';

% Fill in extra args
nArgs = 6;
for ii=nArgs:2:nargin
    optname = varargin{ii-nArgs+1};
    optval = varargin{ii-nArgs+2};
    elike.(optname) = optval;
end

elike.solverOpts = struct();

end

