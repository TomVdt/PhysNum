CXX = clang++

# Executable name (target)
EXE = ex3
CXXFLAGS = -std=c++17 -Wall -Wextra -Wuninitialized -O2
BIN = bin
SRC = src

all: dir $(BIN)/$(EXE)

$(BIN)/$(EXE): $(SRC)/Ex3_2024.cpp
	@$(CXX) -o $@ $^ $(CXXFLAGS) $(DEFINES)

dir:
	@mkdir -p data

clean:
	@rm $(BIN)/$(EXE) || true
	@rm -r data/*.out || true
