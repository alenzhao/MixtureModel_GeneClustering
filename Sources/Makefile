# A Makefile for GITTAR
objects = Gibbs.o mysub.o mathsub.o CmdLine.o

gittar: $(objects)
	g++ -o gittar $(objects) 

Gibbs.o mysub.o mathsub.o: mysub.h
	g++ -O3 -c Gibbs.cpp mysub.cpp mathsub.cpp
CmdLine.o: CmdLine.h
	g++ -O3 -c CmdLine.cpp

clean:
	rm -f $(objects)
	rm -f gittar gittar_without_seq

