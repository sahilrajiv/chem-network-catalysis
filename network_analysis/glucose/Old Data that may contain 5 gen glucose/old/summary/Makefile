SHELL=/bin/bash
summary/all.tex: gen
	echo '\summaryRule{Cannizarro 2, glucose (oxidation)}{\summaryRuleSpan{out/001_r_8_10300000}{out/001\_r\_8\_10300000}}{\summaryRuleSpan{out/002_r_8_11300100}{out/002\_r\_8\_11300100}}{\summaryRuleCombined{out/003_r_8_combined}{out/003\_r\_8\_combined}}{out/004_r_8_constraints.tex}' >> summary/all.tex
	echo '\summaryChapter{Alkaline Glucose Degradation}' >> summary/all.tex
summary.pdf: summary/all.tex
	/usr/local/bin/../bin/mod_post --mode compileSummary
include summary/MakefileImpl
include summary/MakefileDep
include summary/MakefileAllTarget
include summary/MakefileClean
