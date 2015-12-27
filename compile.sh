bibtex main
pdflatex -synctex=1 -interaction=nonstopmode --src-specials main
pdflatex -synctex=1 -interaction=nonstopmode --src-specials main
#for file in main.aux main.bbl main.blg main.log main.lof main.out main.toc
#do
#  mv $file tmp/
#done
