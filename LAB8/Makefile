
CXX = g++
CXXFLAGS = -std=c++17 -I./include -Wall -Wextra -pedantic -O2
LDFLAGS = 

SRC_DIR = src
OBJ_DIR = build/obj
BIN_DIR = build

SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
EXECUTABLE = $(BIN_DIR)/population_simulation


all: directories $(EXECUTABLE)


directories:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

# Компиляция исполняемого файла
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo "Сборка завершена. Исполняемый файл: $@"

# Компиляция объектных файлов
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Очистка
clean:
	rm -rf build/*
	@echo "Очистка завершена."

# Пересборка
rebuild: clean all

# Запуск
run: all
	@echo "Запуск программы..."
	@./$(EXECUTABLE)

# Отладочная сборка
debug: CXXFLAGS += -g -DDEBUG
debug: clean all
	@echo "Отладочная сборка завершена."

# Статический анализ
check:
	cppcheck --enable=all --std=c++17 --suppress=missingIncludeSystem -I./include ./src

# Форматирование кода
format:
	clang-format -i $(SRC_DIR)/*.cpp include/*.h

.PHONY: all clean rebuild run debug check format directories