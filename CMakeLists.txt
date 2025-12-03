# ========================================================
# ПРОЕКТ: Trigger
# БИБЛИОТЕКИ: Eigen, Boost, GSL, matplotlib-cpp
# ========================================================

cmake_minimum_required(VERSION 3.16)  # Повышаем для лучшей поддержки Python3
project(Trigger VERSION 1.0.0 LANGUAGES CXX)

# ========================================================
# НАСТРОЙКИ КОМПИЛЯЦИИ
# ========================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Автоматическое определение типа сборки
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ========================================================
# 1. ПРОВЕРКА И НАСТРОЙКА PYTHON (ВЫПОЛНЯЕТСЯ ПЕРВОЙ!)
# ========================================================
message(STATUS "=== НАСТРОЙКА PYTHON ===")

# Явно задаем нужную версию Python (3.12 в вашем случае)
set(Python3_FIND_VERSION 3.12)
set(Python3_FIND_VERSION_EXACT TRUE)

# Ищем Python с компонентами
find_package(Python3 3.12 REQUIRED COMPONENTS Interpreter Development NumPy)

if(Python3_FOUND)
    message(STATUS "Python3 найден: ${Python3_VERSION}")
    message(STATUS "  Интерпретатор: ${Python3_EXECUTABLE}")
    message(STATUS "  Включения: ${Python3_INCLUDE_DIRS}")
    message(STATUS "  Библиотеки: ${Python3_LIBRARIES}")
    message(STATUS "  NumPy: ${Python3_NumPy_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "Python 3.12 не найден! Установите: sudo apt install python3.12-dev python3-numpy")
endif()

# Резервный поиск NumPy если CMake не нашел
if(NOT Python3_NumPy_INCLUDE_DIR)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_VARIABLE NUMPY_ERROR
    )
    if(NUMPY_ERROR)
        message(FATAL_ERROR "NumPy не найден! Установите: pip3 install numpy")
    endif()
    set(Python3_NumPy_INCLUDE_DIR ${NUMPY_INCLUDE_DIR})
    message(STATUS "  NumPy (ручной поиск): ${Python3_NumPy_INCLUDE_DIR}")
endif()

# ========================================================
# 2. ПОИСК ОСТАЛЬНЫХ БИБЛИОТЕК
# ========================================================
message(STATUS "=== ПОИСК БИБЛИОТЕК ===")

# Eigen3 (заголовочная библиотека)
find_package(Eigen3 3.3 REQUIRED)
message(STATUS "Eigen3: ${EIGEN3_INCLUDE_DIR}")

# Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.65 REQUIRED COMPONENTS filesystem system program_options)
message(STATUS "Boost: ${Boost_INCLUDE_DIRS}")

# GSL
find_package(GSL REQUIRED)
message(STATUS "GSL: ${GSL_INCLUDE_DIRS}")

# ========================================================
# 3. MATPLOTLIB-CPP КОНФИГУРАЦИЯ
# ========================================================
message(STATUS "=== НАСТРОЙКА MATPLOTLIB-CPP ===")

# Проверяем наличие matplotlibcpp.h
if(NOT EXISTS "/usr/local/include/matplotlibcpp.h")
    message(WARNING "matplotlibcpp.h не найден в /usr/local/include/")
    message(STATUS "Скачиваю matplotlib-cpp...")
    
    # Скачиваем и копируем локально
    file(DOWNLOAD 
        https://raw.githubusercontent.com/lava/matplotlib-cpp/master/matplotlibcpp.h
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/matplotlibcpp.h
    )
    set(MATPLOTLIB_CPP_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/libs)
else()
    set(MATPLOTLIB_CPP_INCLUDE /usr/local/include)
endif()

message(STATUS "matplotlibcpp.h: ${MATPLOTLIB_CPP_INCLUDE}")

# Создаем интерфейсную цель для matplotlib-cpp
add_library(matplotlib_cpp INTERFACE)
target_include_directories(matplotlib_cpp INTERFACE
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIR}
    ${MATPLOTLIB_CPP_INCLUDE}
)

# ========================================================
# 4. ИСПОЛНЯЕМЫЙ ФАЙЛ
# ========================================================
message(STATUS "=== СОЗДАНИЕ ИСПОЛНЯЕМОГО ФАЙЛА ===")

# Файлы проекта
file(GLOB SOURCE_FILES "src/*.cpp")
file(GLOB HEADER_FILES "include/*.h")

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

# Включения для нашего проекта
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ========================================================
# 5. ПОДКЛЮЧЕНИЕ ВСЕХ БИБЛИОТЕК
# ========================================================
message(STATUS "=== ПОДКЛЮЧЕНИЕ БИБЛИОТЕК ===")

target_link_libraries(${PROJECT_NAME}
    # Eigen3 (заголовочная)
    Eigen3::Eigen
    
    # Boost
    Boost::filesystem
    Boost::system
    Boost::program_options
    
    # GSL
    GSL::gsl
    GSL::gslcblas
    
    # matplotlib-cpp
    matplotlib_cpp
    
    # Python3
    Python3::Python
    
    # Системные библиотеки
    pthread     # Многопоточность (требуется Python)
    dl          # Динамическая загрузка (требуется Python)
    m           # Математическая библиотека
)

# ========================================================
# 6. НАСТРОЙКИ КОМПИЛЯЦИИ И СБОРКИ
# ========================================================
message(STATUS "=== НАСТРОЙКИ КОМПИЛЯЦИИ ===")

# Флаги компилятора
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter
    )
    
    # Флаги для отладки/релиза
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${PROJECT_NAME} PRIVATE -g -O0)
        target_compile_definitions(${PROJECT_NAME} PRIVATE _DEBUG)
    else()
        target_compile_options(${PROJECT_NAME} PRIVATE -O3 -DNDEBUG)
    endif()
    
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

# ========================================================
# 7. ПОЛЬЗОВАТЕЛЬСКИЕ ЦЕЛИ И УТИЛИТЫ
# ========================================================

# Цель для очистки
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
    COMMENT "Очистка директории сборки"
)

# Вывод информации о проекте
add_custom_target(info
    COMMAND ${CMAKE_COMMAND} -E echo "=== ИНФОРМАЦИЯ О ПРОЕКТЕ ==="
    COMMAND ${CMAKE_COMMAND} -E echo "Проект: ${PROJECT_NAME}"
    COMMAND ${CMAKE_COMMAND} -E echo "Тип сборки: ${CMAKE_BUILD_TYPE}"
    COMMAND ${CMAKE_COMMAND} -E echo "Компилятор: ${CMAKE_CXX_COMPILER}"
    COMMAND ${CMAKE_COMMAND} -E echo "Стандарт C++: ${CMAKE_CXX_STANDARD}"
    COMMENT "Вывод информации о проекте"
)

# ========================================================
# 8. УСТАНОВКА
# ========================================================
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# ========================================================
# 9. ФИНАЛЬНОЕ СООБЩЕНИЕ
# ========================================================
message(STATUS "")
message(STATUS "========================================")
message(STATUS "КОНФИГУРАЦИЯ ЗАВЕРШЕНА УСПЕШНО")
message(STATUS "Проект: ${PROJECT_NAME}")
message(STATUS "Тип сборки: ${CMAKE_BUILD_TYPE}")
message(STATUS "Исполняемый файл: ${PROJECT_NAME}")
message(STATUS "Для сборки выполните: make")
message(STATUS "Для очистки: make clean-all")
message(STATUS "Для информации: make info")
message(STATUS "========================================")