CXX = clang++

# Executable name (target)
EXE = ex4
CXXFLAGS = -std=c++17 -Wall -Wextra -Wuninitialized -O2
BIN = bin
SRC = src

all: dir $(BIN)/$(EXE)

$(BIN)/$(EXE): $(SRC)/Exercice4_2024_students.cpp
	@$(CXX) -o $@ $^ $(CXXFLAGS) $(DEFINES)

dir:
	@mkdir -p data

clean:
	@rm $(BIN)/$(EXE) || true
	@rm -r data/*.out || true
