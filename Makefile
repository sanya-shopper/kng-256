PDF = main.pdf

all: $(PDF)

# Regenerate the last-updated stamp (Pacific time) before every build.
stamp:
	printf '\\newcommand{\\lastupdated}{%s}\n' \
	  "$$(TZ=America/Los_Angeles date '+%B %-d, %Y, %-I:%M %p Pacific')" \
	  > lastupdated.tex

$(PDF): main.tex preamble.tex backmatter.tex sections/*.tex bib/*.bib tools/invert_sigma.py stamp
	latexmk -pdf -interaction=nonstopmode main.tex

watch:
	latexmk -pdf -pvc -interaction=nonstopmode main.tex

clean:
	latexmk -C
	rm -f main.brf main.bbl lastupdated.tex

.PHONY: all watch clean stamp
