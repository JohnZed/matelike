function [cintL, cintU] = elConfRegion(elres, whichVars, confVals, doPlot, plotOpt)
% Computes an EL likelihood-ratio confidence interval and possibly plots it.
% Currently only supports univariate confidence intervals.
%
% Primary Usage:
%
%    "elConfRegion(elres, [whichVars], [confVals], [doPlot], [plotOpt])"
%
%        Arguments in [ ] are optional.
%
%
% Alternative Usage:
%
%    "elConfRegion(elres, whichVars, confVals, 'precise')"
%
%        Slightly more accurate, but much slower, optimization-based means of
%        finding confidence regions.
%
% Inputs:
%
%    "elres" is a result structure previously returned by elSolve(...)
%
%    "whichVars" specifies the indices of theta for which confidence regions
%        are desired. Leave as '[]' or omit to compute CRs for all indices.
%
%    "confVals" specifies one or more confidence levels (e.g. [0.95,0.99]) to
%        consider. Defaults to 0.95
%
%    "doPlot" specifies that the profile p-value should be plotted
%
%    "plotOpt" is an optional structure with fields specifying plotting
%        options. Options include "plotOpt.numPoints" (number of points to sample),
%        "plotOpt.varnames" (a cell array of variable names), "plotOpt.xciLabel"
%        (whether to label the confidence interval on the x-axis), and
%        "plotOpt.horizCR" (whether to draw a horizontal line for the
%        confidence interval).
%
  nv = length(elres.theta);
  
  if ~exist('whichVars') || isempty(whichVars); whichVars = 1:nv; end;
  if ~exist('doPlot') || isempty(doPlot); doPlot = false; end;
  if ~exist('plotOpt') || isempty(plotOpt); plotOpt = struct(); end;
  if ~exist('confVals') || isempty(confVals); confVals = 0.95; end;
  
  % Precise CIs use optimization to find CI crossing instead of grid approach
  if ischar(doPlot) && strcmpi(doPlot,'precise')
    preciseCI = true;
  else
    preciseCI = false;
  end
  
  % Plotting options structure
  npt = fieldopt(plotOpt,'numPoints',10);
  varnames = fieldopt(plotOpt,'varnames','');
  smoothMeth = fieldopt(plotOpt,'smoothMeth','pchip');
  xciLabel = fieldopt(plotOpt,'xciLabel',true);
  horizCR = fieldopt(plotOpt,'horizCR',false);
  rotateLabels = fieldopt(plotOpt,'rotateLabels',false);
  plotPoints = fieldopt(plotOpt,'plotPoints',false);
  labelCI = fieldopt(plotOpt,'labelCI',true);
  adaptiveGrid = fieldopt(plotOpt,'adaptiveGrid',true);
  
  nchangev = length(whichVars);
  alpha = min(confVals, 1 - confVals);
  
  % Store the optimized values
  theta0 = elres.theta;
  [ll0,lr0] = elValues(elres,theta0,[]);
  pv0 = 1;

  % Decide what region to search for the CR
  % Look in a region twice as wide as asymptotic stderr implies
  stderrs = elModelSumm(elres,[],true);
  searchSize = -2 * norminv(min(alpha)) * stderrs;
  minTheta = theta0 - searchSize;
  maxTheta = theta0 + searchSize;
  
  cintL = zeros(nchangev, length(confVals));
  cintU = zeros(nchangev, length(confVals));  
  
  % Layout for multiple graphs
  nrow = floor(sqrt(nchangev));
  ncol = ceil(nchangev / nrow);
  
  % Compute CI for each variable
  for iv=1:nchangev
    vidx = whichVars(iv);
    lb = minTheta(vidx);
    ub = maxTheta(vidx);
    assert(lb < ub);
    
    if preciseCI
      % If we're using the precise method, skip all of the
      % grid and plotting stuff that follows
      for ic=1:numel(confVals)
        thisAlpha = alpha(ic);
        [thisL,thisU] = ciSearch(vidx, lb, ub, alpha, true);
        cintL(iv,ic) = thisL;
        cintU(iv,ic) = thisU;
      end
      
      continue;
    end
    
    if doPlot; subplot(nrow,ncol,iv); end;
    fixv = logical(zeros(nv,1));
    fixv(vidx) = true;
    
    % First use a coarse grid
    coarsePt = npt;
    xlist1 = repmat(theta0,1,coarsePt);
    xlist1(vidx,:) = linspace(lb,ub,coarsePt);
    [prof_ll1,prof_lr1,prof_pv1] = elValues(elres, xlist1, fixv, 1);
    
    
    if adaptiveGrid      
      % Experimental: use a finer grid in the region of interest
      
      % Find the region of interest (near CR boundaries)
      pvmin = min(alpha)/3; pvmax = max(alpha)*3;
      regIdx = find(prof_pv1 >= pvmin & prof_pv1 <= pvmax);
      regPts = xlist1(vidx,regIdx);
    
      % Compute the added grid points
      oldstep = (ub - lb) / (coarsePt - 1);
      newstep = oldstep / 3;
      finePts = sort([regPts - newstep, regPts + newstep]);
      numFine = length(finePts);

      xlist2 = repmat(theta0,1,numFine);
      xlist2(vidx,:) = finePts;
      [prof_ll2,prof_lr2,prof_pv2] = elValues(elres,xlist2,fixv, 1);
    else
      prof_ll2 = []; prof_lr2= []; prof_pv2 = [];
      xlist2 = zeros(length(theta0),0);
    end

    % Combine the two grids and the original estimate
    xlist = [xlist1, xlist2, theta0];
    pvlist = [prof_pv1, prof_pv2, pv0];
    lrlist = [prof_lr1, prof_lr2, lr0];
    % Keep only sorted, unique values
    [tmp,xidx] = union(xlist(vidx,:),[]);
    xlist = xlist(:,xidx);
    pvlist = pvlist(:,xidx);
    lrlist = lrlist(:,xidx);    
    thisxlist = xlist(vidx,:);
    
    % Interpolate on a finer grid to get a smoother plot
    fineX = sort([linspace(lb,ub,npt*10), theta0(vidx)]);
    fineY = interp1(thisxlist, pvlist, fineX, smoothMeth);
    
    for ic=1:numel(confVals)
      [thisL,thisU] = ciSearch(vidx, lb, ub, alpha(ic), false, thisxlist, lrlist);
      cintL(iv,ic) = thisL;
      cintU(iv,ic) = thisU;
    end
    
    if doPlot
      % Plot the p-value profiles and mark confidence regions
      if plotPoints
        plot(fineX, fineY, 'b-', thisxlist, pvlist, 'b.');
      else
        plot(fineX, fineY, 'b-');
      end
      for ic=1:numel(confVals)
        % Choose a color different from that used for main lines
        colorlist = get(gca,'ColorOrder');
        thiscolidx = mod(ic,size(colorlist,1)-1)+2;
        thiscol = colorlist(thiscolidx,:);
        
        if horizCR
          line([lb,ub], [alpha(ic),alpha(ic)], 'LineStyle','--','Color',thiscol);
        end
        
        % Draw vertical lines for the region
        line([cintL(iv,ic),cintL(iv,ic)],[0,1],...
             'LineStyle','--','Color',thiscol);
        line([cintU(iv,ic),cintU(iv,ic)],[0,1],...
             'LineStyle','--','Color',thiscol);
        % Optionally label the lines
        if labelCI
          vertPos = 0.9;
          tstr = sprintf('%4.3g%% CI', 100*(1-alpha(ic)));
          txtOpts = {'Color',thiscol,'Rotation',-90};
          hText1 = text(cintL(iv,ic), vertPos, tstr, txtOpts{:},'VerticalAlignment','top');
          hText2 = text(cintU(iv,ic), vertPos, tstr, txtOpts{:},'VerticalAlignment','bottom');
        end
      end
      xlim([lb,ub]);
      ylim([0,1]);
      if isempty(varnames)
        title(sprintf('var %d',vidx));
      else
        title(varnames{iv});
      end
      ylabel('P-value');
      
      % Mark the CI points on both axes
      ylnow = get(gca,'ylim');
      numticks = 4;
      baseticky = linspace(max(alpha)+0.1,1, numticks);
      tvy = sort([baseticky, alpha]);
      tnamey = num2str(tvy','%4.2f');
      set(gca,'YTick',tvy,'YTickLabel',tnamey);
      
      if xciLabel
        % Manually add xticks for the boundaries of the confidence region
        xlnow = get(gca,'xlim');
        xtnames = [1 - alpha, 1 - alpha] * 100;
        basetickx = linspace(minTheta(iv),maxTheta(iv),numticks);
        newtickx = [theta0(iv), cintL(iv,:), cintU(iv,:)];
        % Drop too-close values
        xthresh = 0.05*(xlnow(2) - xlnow(1));
        [XT1,XT2] = meshgrid(basetickx,newtickx);
        badtick = any(abs(XT1 - XT2) < xthresh,1);
        tvx = sort([basetickx(~badtick),newtickx]);
        
        set(gca,'XTick',tvx,'XTickLabel','');  % Remove old labels
        
        % Rotation trick inspired by xticklabel_rotate script by BFG Katz
        xtlbl = num2str(tvx','%5.3f');
        xtPos = get(get(gca,'XLabel'),'Position');
        hText = text(tvx, repmat(xtPos(2),size(tvx,1),size(tvx,2)), xtlbl);
        
        if rotateLabels
          set(hText,'Rotation',-90,'HorizontalAlignment','left');
        end
      end
    end
  end
  
  
function [thisL,thisU] = ciSearch(vidx, lb, ub, ciAlpha, precise, xlist, lrlist)
  fzopts = optimset('fzero');
  fzopts.TolX = 1e-12;
  % Note: precision of solution is actually much worse than TolX, as the
  % underlying optimization in elValues adds error
  
  % First search below thetahat, then above it to get the
  % two boundaries of the confidence intervals
  
  if precise
    % Precise method recomputes profile LR each time
    target = chi2inv(1 - ciAlpha, 1);
    thisFun = @(x) preciseCIobj(x, vidx, target);
  else
    % Interpolate through list of p-values (faster, but less accurate)
    target = chi2inv(1 - ciAlpha, 1);
    thisFun = @(x) interp1(xlist,lrlist,x,smoothMeth) - target;
  end
  
  % Lower bound
  thisReg = [lb,theta0(vidx)];
  [thisL,fval,flag1] = fzero(thisFun, thisReg, fzopts);
  
  % Upper bound
  thisReg = [theta0(vidx), ub];
  [thisU,fval,flag2] = fzero(thisFun, thisReg, fzopts);
  
  if flag1 <= 0 || flag2 <= 0
    fprintf('Failed to solve inner EL problem for CI: idx=%d, alpha=%g\n',...
            vidx, ciAlpha);
    thisL = NaN;
    thisU = NaN;
  end
end

function f = preciseCIobj(x, vidx, target)
% Computes the difference between this LR stat and the desired one
  if x == theta0(vidx)
    % Optimization: just return the value if we know it    
    f = lr0 - target;
    return;
  end
    
  thisList = theta0;
  thisList(vidx) = x;
  [llOut,lrOut] = elValues(elres, thisList, vidx);
  f = lrOut - target;
end

end