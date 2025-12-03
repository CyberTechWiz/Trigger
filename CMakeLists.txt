# Минимальная версия CMake
cmake_minimum_required(VERSION 3.10)

# Имя проекта
project(Trigger VERSION 1.0.0 LANGUAGES CXX)

# Используем стандарт C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Требуем именно этот стандарт
set(CMAKE_CXX_EXTENSIONS OFF)        # Отключаем расширения компилятора

# === ПОИСК И НАСТРОЙКА БИБЛИОТЕК ===

# 1. EIGEN (заголовочная библиотека - только .h файлы)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
# NO_MODULE - ищем через конфиг файл, а не через FindEigen.cmake
# REQUIRED - прервёмся если не найдём

# 2. BOOST (ищем конкретные компоненты)
find_package(Boost 1.65 REQUIRED COMPONENTS filesystem system program_options)
# COMPONENTS - какие части Boost нужны:
# - filesystem: работа с файлами
# - system: системные ошибки
# - program_options: парсинг аргументов командной строки

# 3. GSL (GNU Scientific Library)
find_package(GSL REQUIRED)
# Найдёт libgsl и libgslcblas

# 4. Python для matplotlib-cpp
find_package(PythonLibs 3.6 REQUIRED)
# Найдёт Python.h и библиотеки Python

# 5. matplotlib-cpp (просто заголовочный файл)
# Создаём "интерфейсную" цель - она только указывает пути
add_library(matplotlibcpp INTERFACE)
target_include_directories(matplotlibcpp INTERFACE
    ${PYTHON_INCLUDE_DIRS}            # Путь к Python.h
    "/usr/local/include"              # Где лежит matplotlibcpp.h
    # Альтернативно можно положить matplotlibcpp.h в папку проекта
    # "${CMAKE_CURRENT_SOURCE_DIR}/libs"
)

# === СОЗДАНИЕ ИСПОЛНЯЕМОГО ФАЙЛА ===

# Все исходные файлы программы
set(SOURCES
    src/main.cpp
    # src/math_utils.cpp
    # src/plotter.cpp
)

# Все заголовочные файлы
set(HEADERS
    include/math_utils.h
    # include/plotter.h
)

# Создаём исполняемый файл
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Указываем, где искать заголовочные файлы ВАШЕГО проекта
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# === ПОДКЛЮЧЕНИЕ БИБЛИОТЕК К ПРОГРАММЕ ===

target_link_libraries(${PROJECT_NAME}
    # Eigen (только заголовки, но указываем для порядка)
    Eigen3::Eigen
    
    # Boost компоненты
    Boost::filesystem
    Boost::system
    Boost::program_options
    
    # GSL библиотеки
    GSL::gsl      # Основная библиотека
    GSL::gslcblas # BLAS совместимость
    
    # matplotlib-cpp
    matplotlibcpp
    
    # Python библиотеки
    ${PYTHON_LIBRARIES}
    
    # Математическая библиотека (часто нужна)
    m
)

# === ДОПОЛНИТЕЛЬНЫЕ НАСТРОЙКИ ===

# Включить все предупреждения компилятора
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE
        -Wall        # Все стандартные предупреждения
        -Wextra      # Дополнительные предупреждения
        -Wpedantic   # Строгое соответствие стандарту
        -Werror      # Превратить предупреждения в ошибки (опционально)
    )
endif()

# Оптимизация для релиза
set(CMAKE_BUILD_TYPE Release CACHE STRING "Тип сборки" FORCE)

# Установка программы (опционально)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)