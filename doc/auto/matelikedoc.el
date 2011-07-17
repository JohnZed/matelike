(TeX-add-style-hook "matelikedoc"
 (lambda ()
    (LaTeX-add-bibliographies
     "elpapers")
    (LaTeX-add-labels
     "sec:introduction"
     "sec:simple-example"
     "eq:momExample"
     "sec:empirical-likelihood"
     "sec:basic-el-model"
     "eq:elmax"
     "eq:crmax"
     "sec:el-computation"
     "sec:installing-melike"
     "sec:basic-melike"
     "sec:sample-results"
     "sec:linear-iv"
     "tab:dpd-res"
     "sec:poisson-iv-model"
     "tab:cig-res"
     "changes"
     "sec:references")
    (TeX-add-symbols
     '("linesect" 1)
     '("R" 1)
     '("DAdxy" 3)
     '("Ddxy" 2)
     '("Dydx" 2)
     '("Dydxx" 2)
     '("Ddxx" 1)
     '("Ddx" 1)
     '("oneover" 1)
     '("Exb" 1)
     "melike"
     "sumin"
     "half")
    (TeX-run-style-hooks
     "natbib"
     "url"
     "listings"
     "amsfonts"
     "amsmath"
     "setspace"
     "latex2e"
     "art10"
     "article"
     "mellisting")))

