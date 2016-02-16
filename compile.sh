echo "#################################"
echo "#################################"
echo "#################################"
cp /Users/lacerda/Documents/Papers/papers.bib bibliografia.bib
echo "#################################"
echo "#################################"
echo "#################################"
pdflatex -synctex=1 -interaction=nonstopmode --src-specials main
echo "#################################"
echo "#################################"
echo "#################################"
bibtex main
echo "#################################"
echo "#################################"
echo "#################################"
pdflatex -synctex=1 -interaction=nonstopmode --src-specials main
echo "#################################"
echo "#################################"
echo "#################################"
pdflatex -synctex=1 -interaction=nonstopmode --src-specials main
#for file in main.aux main.bbl main.blg main.log main.lof main.out main.toc
#do
#  rm -f $file
#done
