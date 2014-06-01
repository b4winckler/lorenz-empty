all: hyperbolic.pdf crenorm.pdf

empty.txt: EmptyRenormalization.hs
	runhaskell EmptyRenormalization.hs 1000 > empty.txt

%.pdf: empty.txt gengraphs.R
	Rscript gengraphs.R
