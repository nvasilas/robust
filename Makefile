MAIN = main
CHAPTERS = chapters
BIB = bible

LATEX=rubber -f --pdf --synctex -s
XELATEX=rubber -f --pdf --synctex -s --module xelatex
BUILD=$(LATEX)

DEPENDS = $(wildcard *.tex)
ifneq ($(CHAPTERS),)
	DEPENDS += $(wildcard $(CHAPTERS)/*.tex)
endif
ifneq ($(BIB),)
	DEPENDS += $(BIB).bib
endif

all: $(MAIN).pdf

%.pdf: %.tex $(DEPENDS)
	$(BUILD) $<
	rubber-info --check $<

clean:
	rubber --clean $(DEPENDS)
	rm -f $(MAIN).pdf

.PHONY : all clean
