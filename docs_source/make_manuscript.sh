#make latexpdf
make latex
cd _build/latex/
sed -i 's/sphinxhowto/article/g' bionumpymanuscript.tex
tectonic bionumpymanuscript.tex
gio open bionumpymanuscript.pdf || true  # ignore error if gio does not exist
