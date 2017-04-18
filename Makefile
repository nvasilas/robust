MAIN = main
CHAPTERS = chapters
COVER = cover
BIB = bible
PREAMBLE = preamble
MACROS = macros

LATEX=rubber -f --pdf --synctex -s
XELATEX=rubber -f --pdf --synctex -s --module xelatex
BUILD=$(LATEX)

DEPENDS = $(MAIN).tex
ifneq ($(CHAPTERS),)
	DEPENDS += $(wildcard $(CHAPTERS)/*.tex)
endif
ifneq ($(COVER),)
	DEPENDS += $(COVER).tex
endif
ifneq ($(BIB),)
	DEPENDS += $(BIB).bib
endif
ifneq ($(PREAMBLE),)
	DEPENDS += $(PREAMBLE).tex
endif
ifneq ($(MACROS),)
	DEPENDS += $(MACROS).tex
endif

all: $(MAIN).pdf

%.pdf: %.tex $(DEPENDS)
	$(BUILD) $<
	rubber-info --check $<

clean:
	rubber --clean $(DEPENDS)
	rm -f $(MAIN).pdf

.PHONY : all clean
