ALL:	align consensus
align:	align.cc
	g++ -o align align.cc 
consensus: consensus.cc
	g++ -o consensus consensus.cc
clean:	
	rm align consensus coffee* *~ 
