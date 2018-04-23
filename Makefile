C = g++
FLAGS = -Wall -O3 -fopenmp

all: biosig.cpp
	$(C) $(FLAGS) biosig.cpp -o biosig