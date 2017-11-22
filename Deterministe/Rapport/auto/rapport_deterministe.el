(TeX-add-style-hook
 "rapport_deterministe"
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
    "output_1_neg1_5_1"
    "output_1_1_5_1"
    "output_1_05_5_2"
    "output_two_steps_1_1_1_2"
    "output_two_steps_1_05_1_2"
    "output_two_steps_1_05_3_2"
    "output_limit_sigma_a_zero"
    "loop_nb_pts_mu_delta"
    "loop_nb_pts_mu_cste"
    "nb_segs_diff_L2_delta"
    "nb_segs_diff_L2_cste"
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
    '("questions" 2)
    '("question" 2)
    '("norm" 1)
    '("dx" 1)
    "tphi"
    "intsigma"
    "F"
    "Qt"
    "Phit")
   (LaTeX-add-labels
    "eq:principal"
    "eq:pasdiff"
    "eq:recu"
    "fig:exple_a"
    "fig:exple_b"
    "fig:exple_c"
    "fig:ts_1"
    "fig:ts_05"
    "fig:ts_3"
    "eq:eps"
    "fig:limit_sigma_zero"
    "tab:DSA"
    "fig:lm_d"
    "fig:lm_c"
    "fig:Nx_diff_d"
    "fig:Nx_diff_c"))
 :latex)

