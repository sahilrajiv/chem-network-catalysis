SHELL=/bin/bash
summary/all.tex: gen
	echo '\summaryChapter{Formose Reaction}' >> summary/all.tex
summary.pdf: summary/all.tex
	/usr/local/bin/../bin/mod_post --mode compileSummary
include summary/MakefileImpl
include summary/MakefileDep
include summary/MakefileAllTarget
include summary/MakefileClean
