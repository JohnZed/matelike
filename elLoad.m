function elLoad(solver, forceReload)
% Loads ipopt, matelike, and intlab libraries
%
% Usage:
%
%    "elLoad"
% or:
%     "elLoad(solver, forceReload)"
%
% Inputs:
%
%    "solver" is an optional string specifying the nonlinear solver to use.
%
%    Valid solver options are:
%
%        'zipsolver' -- matElike's default, built-in solver
%        'ipopt' -- the Ipopt solver, which require the Ipopt MEX interface and
%                library to be installed
%        'fmincon' -- the fmincon solver from Matlab's optimization
%                toolbox. This is only available with optimization toolbox
%                version 4.0 (Matlab R2008a) or greater, as the older versions
%                lack support for analytic Hessians with inequality constraints.
%
%     In general, Ipopt is the most robust solver, but 'zipsolver' works well
%     for most problems as long as they're reasonably well scaled. 'fmincon'
%     appears to be less reliable, but slightly faster than zipsolver.
% 
% IMPORTANT: You will need to edit BASEPATH below to point to the directory
% containing the 'matelike'  and 'zipsolver' folders on your system.
%

BASEPATH = '~/matlab/melroot';

% You don't have to modify the following if you use the default install
INTLABPATH =    [BASEPATH '/intlab/'];
ZIPSOLVERPATH = [BASEPATH '/zipsolver/'];
ELIKEPATH =     [BASEPATH '/matelike/'];
IPOPTPATH = [BASEPATH '/ipopt/'];

% Set to false if you don't have Intlab, e.g. for licensing reasons
useIntlab = true;

elVersion = '0.1.4';

% Set the solver to use
global elDefaultSolver;
if (nargin >= 1)
  elDefaultSolver = solver;
elseif ~isempty(elDefaultSolver)
  % Keep the current default
else
  % Default to zipsolver
  elDefaultSolver = 'zipsolver';
end

if verLessThan('matlab','7.4')
  fprintf('Warning: matElike not tested with your version of Matlab\n');
end

global elLoaded;
if (isempty(elLoaded) || (nargin >=2 && forceReload))
  addpath([IPOPTPATH 'matlab']);
  addpath(ZIPSOLVERPATH);
  addpath(ELIKEPATH);
  if (useIntlab)
    loadIntlabAD;
  end
  elLoaded = 1;
end


function loadIntlabAD()
% Only load the automatic differentiation component from Intlab
% This is based on startintlab.m from the Intlab 5.4 distribution

disp('Loading Intlab Automatic differentiation');

% Banner suggestd by Siegfried Rump as condition for including Intlab in
% matElike distribution (personal email to JohnZ, Tues Mar 11, 2008).
fprintf('%s\n',...
    '=========================================================================  ',... 
     '***** Uses INTLAB - INTerval LABoratory      www.ti3.tu-harburg.de/ ***** ',... 
     '*****   The Matlab toolbox for Reliable Computing                   ***** ',... 
     '*****   (c) Siegfried M. Rump, Insitute for Reliable Computing      ***** ',... 
     '*****   Hamburg University of Technology, Germany                   ***** ',... 
     '========================================================================= ',... 
     '*****  Free for private and academic use. Commercial use or any use ***** ',... 
     '*****  in conjunction with a commercial program requiring this      ***** ',... 
     '*****  program or part of it to function properly is prohibited.    ***** ',... 
     '=========================================================================');

addpath( [ INTLABPATH ]                           , ...
         [ INTLABPATH 'gradient' ]                , ...
         [ INTLABPATH 'hessian' ]                 , ...
         [ INTLABPATH 'intval' ]                 , ...         
         [ INTLABPATH 'slope' ]                   , ...
         [ INTLABPATH 'utility' ]                 , ...
         [ INTLABPATH 'long' ]);

  spparms('autommd',0);
  
  gradientinit;
  hessianinit;

  % Always use sparse hessian but dense gradient
  sparsehessian(0);
  sparsegradient(Inf);

  fprintf('\nUsing %s nonlinear solver\n', elDefaultSolver);
  
  fprintf('\n----  Loaded matElike verson %s package successfully ----\n', elVersion);
end



end