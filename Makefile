CXX ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -fopenmp

SRC := src/main.cpp src/graph.cpp src/simulate.cpp src/influence.cpp src/realizations.cpp src/io.cpp
OBJ := $(SRC:.cpp=.o)

build/sim: $(OBJ)
	mkdir -p build
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf build src/*.o
