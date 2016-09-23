CXX = g++
CXXFLAGS = -march=native -std=c++11 -g -O3 -fopt-info
LDFLAGS = -lm -lfftw3 -fopenmp

EXECUTABLE_NAME=main

# Folders
BIN = ../Bin
OBJ = $(BIN)/Obj

# Files
SOURCE_FILES=\
    main.cpp \
    Ket.cpp \
    Representation.cpp

INPUT_FILES = input.dat
EXECUTABLE_FILES = $(EXECUTABLE_NAME:%=$(BIN)/%)
OBJECT_FILES     = $(SOURCE_FILES:%.cpp=$(OBJ)/%.o)
# ^^^ A more succinct expression for $(OBJECT_FILES), using
#     http://www.gnu.org/software/make/manual/make.html#Substitution-Refs

build: $(EXECUTABLE_FILES) input

clean:
	rm -rf $(BIN)

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean

$(EXECUTABLE_FILES): $(OBJECT_FILES)
	@$(CXX) -o $@ $^ $(LDFLAGS)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
	@echo "Build successful!"

# http://www.gnu.org/software/make/manual/make.html#Static-Pattern
$(OBJECT_FILES): $(OBJ)/%.o: %.cpp common.h
	@echo Compiling $<
	@mkdir -p $(@D)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#index-_0024_0028_0040D_0029
	@$(CXX) -c $(CXXFLAGS) $< -o $@

input: $(BIN)/$(INPUT_FILES)
$(BIN)/$(INPUT_FILES): $(INPUT_FILES)
		@cp $(INPUT_FILES) $(BIN)
		@echo "Copied input"
