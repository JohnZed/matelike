% Converts function help docs to latex

% Functions to includes
funclist = { 'elLoad', 'elSetup', 'elSolve', 'elLinearIV', 'elModelSumm', ...
              'elEval', 'elLRTest', 'elConfRegion', 'elValues' };

% Convert to subsection headings
highlight = { 'Inputs', 'Outputs', 'Output', 'Primary Usage', 'Alternative Usage', ...
              'Details', 'Usage', 'Requirements' };
mytexify = @(x) sprintf('\\\\end{lstlisting}\n{\\\\small\\\\textbf{%s:}}\n\\\\begin{lstlisting}',x);
highlightTex = cellfun(mytexify, highlight,...
                       'UniformOutput', false);
highlightSearch = cellfun(@(x) sprintf('\n *%s:',x), highlight, ...
                       'UniformOutput', false);

%fid = 1; % stdout
fid = fopen('mellisting.tex','w');

for ifn=1:length(funclist)
  fn = funclist{ifn};
  fprintf(fid, '\\linesect{%s}\n', fn);
  htext = help(fn);
  
  % Highlight section titles
  newtext = regexprep(htext, highlightSearch, highlightTex, 'ignorecase');

  % Highlight code sections
  newtext = regexprep(newtext, '"([^"\n]+)"', '"\\texttt{$1}"');
  
  fprintf(fid, '\\begin{lstlisting}\n');
  fprintf(fid, '%s\n', newtext);
  
  fprintf(fid, '\\end{lstlisting}\n');
end

fclose(fid);
