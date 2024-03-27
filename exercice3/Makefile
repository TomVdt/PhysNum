CXX = clang++

# Executable name (target)
EXE = ex3
CXXFLAGS = -std=c++23 -Wall -Wextra -O2
BIN = bin
SRC = src

$(BIN)/$(EXE): $(SRC)/Ex2_2024.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) $(DEFINES)

clean:
	@rm $(BIN)/$(EXE) || true
	@rm -r data/*.out || true

# .PHONY clean