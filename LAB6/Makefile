CXX = g++
CXXFLAGS = -std=c++17 -Wall -Iheaders
TARGET = build/main
SOURCES = src/complex_operations.cpp src/data_exporters.cpp src/signal_filters.cpp src/signal_generator.cpp src/transform_analyzers.cpp src/utilities.cpp main.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -rf build

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run