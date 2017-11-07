(TeX-add-style-hook
 "rapport_monte_carlo"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("babel" "francais") ("fontenc" "T1") ("hyperref" "colorlinks=false") ("geometry" "top=2cm" "bottom=2cm" "left=2cm" "right=2cm")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "inputenc"
    "babel"
    "fontenc"
    "amsmath"
    "amsfonts"
    "amssymb"
    "hyperref"
    "geometry"
    "graphicx"
    "caption"
    "color"
    "float"
    "fancyhdr"
    "mathtools")
   (TeX-add-symbols
    '("question" 2)
    '("norm" 1)
    '("dx" 1)
    "tphi"
    "intsigma"
    "F")
   (LaTeX-add-labels
    "eq:principal"
    "fig:delta_many_mu"
    "fig:sigmacst"
    "fig:cste_many_mu"
    "fig:phi"
    "fig:sigmavar"
    "fig:TP1_diff"
    "fig:nbjumps"
    "fig:vcv"))
 :latex)

