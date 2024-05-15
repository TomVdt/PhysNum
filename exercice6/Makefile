CXX = clang++

# Executable name (target)
EXE = ex6
CXXFLAGS = -std=c++17 -Wall -Wextra -Wuninitialized -O2
BIN = bin
SRC = src

all: dir $(BIN)/$(EXE)

$(BIN)/$(EXE): $(SRC)/Exercice6_2024_student.cpp
	@$(CXX) -o $@ $^ $(CXXFLAGS) $(DEFINES)

dir:
	@mkdir -p data

clean:
	@rm $(BIN)/$(EXE) || true
	@rm -r data || true
