(TeX-add-style-hook
 "presentation"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "17pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8")))
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "geometry"
    "inputenc"
    "multirow"
    "amssymb"
    "amsmath"
    "mathtools"
    "graphicx"
    "epstopdf"
    "caption"
    "tikz")
   (TeX-add-symbols
    '("blockincludetwo" 2)
    '("blockinclude" 1)
    "F"
    "U"
    "err"
    "rb")
   (LaTeX-add-labels
    "fig:exple_a"
    "fig:exple_b"
    "fig:ts_05"
    "fig:ts_3"
    "fig:TP1_diff"
    "fig:nbjumps"
    "tab:DSA"
    "fig:vcv"
    "fig:Nx_diff_d"
    "fig:Nx_diff_c"
    "fig:lm_d"
    "fig:lm_c"
    "fig:limit_sigma_zero")
   (LaTeX-add-environments
    '("blockenumerate" 1)
    '("blockitemize" 1))
   (LaTeX-add-lengths
    "halfwidth"
    "halfheight"))
 :latex)

