CXX=g++  
##CCFLAG= -O3 -Wall -Wextra -Werror --std=c++11
CCFLAG= -O3 -Wall  --std=c++14
LDFLAG= -lfftw3 -lfftw3f -lm 

main.x: main.cc 2dSFG.o class.o response2d.o dft.o # class.o
	$(CXX) main.cc -o main.x $(CCFLAG) 2dSFG.o class.o response2d.o dft.o $(LDFLAG)

2dSFG.o: 2dSFG.cc 
	$(CXX) -c $(CCFLAG) 2dSFG.cc 

dft.o: dft.cc 
	$(CXX) -c $(CCFLAG) dft.cc 

response2d.o: response2d.cc #class.o
	$(CXX) -c $(CCFLAG) response2d.cc #class.o

class.o: class.cc
	$(CXX) -c $(CCFLAG) class.cc
	

.PHONY: clean all

clean:
	rm -f *.o *.x
