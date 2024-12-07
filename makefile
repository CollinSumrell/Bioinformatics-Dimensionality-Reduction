# Makefile for PCA Program

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall

# Directories
INCLUDES = -I /opt/homebrew/include/eigen3 -I$(SRC_DIR)/include  # Include the directory where CSVReader.h is located
SRC_DIR = src/cpp-pca
OBJ_DIR = $(SRC_DIR)/obj
BIN_DIR = bin

# Source files and executable
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)  # Automatically include all .cpp files in src directory
OBJ_FILES = $(SRC_FILES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)  # Map .cpp files to .o files
EXEC = $(BIN_DIR)/pca_program

# Targets
all: $(EXEC) run  # Add run as a dependency of the 'all' target

$(EXEC): $(OBJ_FILES)
	@echo "Linking $(OBJ_FILES) to create $(EXEC)"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXEC) $(OBJ_FILES)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $< into $@"
	@mkdir -p $(OBJ_DIR)  # Ensure obj directory exists
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@echo "Cleaning up..."
	rm -rf $(OBJ_DIR)/*.o $(BIN_DIR)/pca_program

run: $(EXEC)
	@echo "Running $(EXEC)..."
	./$(EXEC)

# Phony targets
.PHONY: all clean run