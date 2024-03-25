CXX = clang++

# Executable name (target)
EXE = ex2
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
BIN = bin
SRC = src

$(BIN)/$(EXE): $(SRC)/Ex2_2024_student.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS) $(DEFINES)

clean:
	@rm $(BIN)/$(EXE) || true
	@rm -r data/*.out || true

# .PHONY clean