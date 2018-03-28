C = g++
FLAGS = -Wall -O3

all: biosig.cpp
	$(C) $(FLAGS) biosig.cpp -o biosig