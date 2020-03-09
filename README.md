# DNASequenceAlignment
DNA Sequence alignments with global and local alignment methods.

MAKEFILE CONTENTS:

allalign:
	python alignment.py --mode $(MODE) --input $(INPUT) --gapopen $(GAPOPEN) --gapext $(GAPEXT)

-----------------------------------------------------------------------------------------------------
In order to run the code, type:
	make allaling MODE=<selectedMode> INPUT=<inputFileName> GAPOPEN=<int> GAPEXT=<int>

where the last variable gapExt is necessary only when the <selectedMode> is aglobal or alocal

mode types:
	-global
	-aglobal
	-local
	-alocal

GAPOPEN and GAPEXT: any integer value
-----------------------------------------------------------------------------------------------------
Sample commands to run the code:

	make allalign MODE=global INPUT=sequences.fa GAPOPEN=-5
	make allalign MODE=aglobal INPUT=sequences.fa GAPOPEN=-5 GAPEXT=-2
