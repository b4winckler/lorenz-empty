all: data/hyperbolic.pdf data/crenorm.pdf

data/empty.txt: EmptyRenormalization.hs
	@mkdir -p data
	runhaskell EmptyRenormalization.hs 1000 > data/empty.txt

data/%.pdf: data/empty.txt gengraphs.R
	Rscript gengraphs.R
