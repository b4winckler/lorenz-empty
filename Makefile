all: data/hyperbolic.pdf data/crenorm.pdf

clean:
	rm -rf data

data/empty.txt: EmptyRenormalization.hs
	@mkdir -p data
	runhaskell EmptyRenormalization.hs 1 2 999 > data/empty.txt

data/%.pdf: data/empty.txt gengraphs.R
	Rscript gengraphs.R
