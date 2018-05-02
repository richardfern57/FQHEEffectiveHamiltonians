CC = g++
DEPS = src/Binomial.h src/Factorials.h src/Interactions.h src/utility.h src/Jack.h
OBJ = src/Binomial.o src/Factorials.o src/Interactions.o src/utility.o src/Jack.o src/main.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++0x -c -o $@ $<

Main: $(OBJ)
	$(CC) -std=c++0x -o $@ $^


