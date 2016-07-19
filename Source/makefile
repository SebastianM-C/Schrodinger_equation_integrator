CXX=g++
CXXFLAGS=-march=native -O2 
LDFLAGS = -lm -lfftw3 -fopenmp

EXECUTABLE_NAME=test

# Folders
BIN=../Output/bin
OBJ=$(BIN)/obj

# Files
SOURCE_FILES=\
    main.cpp \
    ket.cpp \
    position.cpp 

EXECUTABLE_FILES = $(EXECUTABLE_NAME:%=$(BIN)/%)
OBJECT_FILES     = $(SOURCE_FILES:%.cpp=$(OBJ)/%.o)
# ^^^ A more succinct expression for $(OBJECT_FILES), using
#     http://www.gnu.org/software/make/manual/make.html#Substitution-Refs

build: $(EXECUTABLE_FILES)

clean:
	rm -r -f $(BIN)
	@# ^^^ I don't recommend suppressing the echoing of the command using @

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean

$(EXECUTABLE_FILES): $(OBJECT_FILES)
	@$(CXX) -o $@ $^ $(LDFLAGS)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
	@cp input.dat $(BIN)
	@echo "Build successful!"

# http://www.gnu.org/software/make/manual/make.html#Static-Pattern
$(OBJECT_FILES): $(OBJ)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#index-_0024_0028_0040D_0029
	@$(CXX) -c $(CXXFLAGS) $< -o $@
