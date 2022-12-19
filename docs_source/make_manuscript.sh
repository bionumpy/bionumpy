make latexpdf
cd _build/latex/
sed -i 's/sphinxhowto/article/g' bionumpymanuscript.tex
pdflatex bionumpymanuscript.tex
gio open bionumpymanuscript.pdf
