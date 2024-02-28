CXX = clang++

# Executable name (target)
EXE = ex1
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
BIN = bin
SRC = src

$(BIN)/ex1: $(SRC)/Exercice1_student.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) $(DEFINES)

clean:
	@rm $(BIN)/ex1
	@rm -r data/*.log || true

# .PHONY clean