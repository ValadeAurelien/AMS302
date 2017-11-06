(TeX-add-style-hook
 "phi"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-add-symbols
    '("includegraphics" ["argument"] 1)
    '("color" ["argument"] 1)
    '("rotatebox" 2)
    "gplgaddtomacro"
    "colorrgb"
    "colorgray")
   (LaTeX-add-lengths
    "gptboxheight"
    "gptboxwidth")
   (LaTeX-add-saveboxes
    "gptboxtext"))
 :latex)

