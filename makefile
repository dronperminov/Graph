COMPILER=g++
FLAGS=-Wall
FILES=main.cpp

all:
	$(COMPILER) $(FLAGS) $(FILES) -o graph

clean:
	rm graph