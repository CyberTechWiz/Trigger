#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <boost/numeric/odeint.hpp>
#include "matplotlibcpp.h"
// Подключаем решатель полиномов из GSL
#include <gsl/gsl_poly.h>
// Библиотеки Boost
#include <boost/numeric/odeint.hpp>
// Библиотека рисования
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace boost::numeric::odeint;

// --- ЧАСТЬ 1: ОПИСАНИЕ СИСТЕМЫ ---
// Это "закон природы" для нашей программы.
// Компьютер будет вызывать этот код тысячи раз, чтобы узнать "куда двигаться дальше".
class BistableSystem
{
    double a; // Параметр сдвига
    double b; // Параметр бифуркации

public:
    // Конструктор: задаем параметры один раз при создании
    BistableSystem(double param_a, double param_b) : a(param_a), b(param_b) {}

    // Оператор (): Это стандарт Boost.
    // x    - вход: текущее значение x (где мы сейчас)
    // dxdt - выход: сюда мы записываем скорость (dx/dt)
    // t    - время (в этом уравнении оно явно не участвует, но нужно для стандарта)
    void operator()(const std::vector<double> &x, std::vector<double> &dxdt, double t)
    {
        // Наше уравнение: dx/dt = a + bx - x^3
        // x[0], потому что вектор может хранить много переменных, но у нас одна.
        dxdt[0] = a + b * x[0] - (x[0] * x[0] * x[0]);
    }
};

// --- ЧАСТЬ 2: НАБЛЮДАТЕЛЬ (ЗАПИСЬ ДАННЫХ) ---
// Boost просто считает. Чтобы сохранить цифры для графика, нужен "Наблюдатель".
// Он "подглядывает" за процессом на каждом шаге и записывает x и t в наши векторы.
struct Observer {
    std::vector<double>& t_vec;
    std::vector<double>& x_vec;

    Observer(std::vector<double>& t, std::vector<double>& x) : t_vec(t), x_vec(x) {}

    void operator()(const std::vector<double> &x, double t) {
        t_vec.push_back(t);
        x_vec.push_back(x[0]);
    }
};

// Структура для хранения информации о найденной точке равновесия
struct EquilibriumPoint {
    double value;       // Координата x
    bool is_stable;     // true = устойчиво, false = неустойчиво
    std::string type;   // Текстовое описание
};

// Функция проверки устойчивости: lambda = b - 3x^2  - производная
EquilibriumPoint analyze_point(double b, double x) {
    double lambda = b - 3.0 * x * x;
    EquilibriumPoint p;
    p.value = x;
    const double eps = 1e-12;  // Маленькое число для сравнения

    if (std::abs(lambda) < eps) {
        p.is_stable = false;  // Полуустойчивые технически не стабильны
        p.type = "полуустойчивое";
    } 

    else if (lambda < 0) {
        p.is_stable = true;
        p.type = "устойчивое";
    } else {
        p.is_stable = false;
        p.type = "неустойчивое";
    }
    return p;
}


// Функция для генерации списка чисел (как в Python numpy.linspace)
// Нам нужно создать диапазон x от -2 до 2, чтобы построить график
std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;
    double step = (end - start) / (num - 1);
    for(int i = 0; i < num; ++i) {
        result.push_back(start + step * i);
    }
    return result;
}

// Анализ устойчивости (как в прошлом уроке)
bool is_stable(double b, double x) {
    double lambda = b - 3.0 * x * x;
    return (lambda < 0);
}

int main()
{
    // 2. Задание значений параметров системы
    double a = 0.0;
    double b = 4.0;
    double x0 = 0.1;
    double t_end = 20;
    BistableSystem system(a, b);

    // 3. Задание начальных условий
    // Начнем с x = 0.1. Система должна "скатиться" в устойчивое состояние
    std::vector<double> x_state = {0.1};

    // Векторы, куда будем писать историю
    std::vector<double> times;
    std::vector<double> x_vals;

    // -------------------------------------------------------------------------------------------------------------------------------
    // 1. Определение точек равновесия и их типа
    // --- ПОДГОТОВКА ДЛЯ GSL ---
    // GSL решает x^3 + c1*x^2 + c2*x + c3 = 0
    // Наше уравнение: x^3 - b*x - a = 0
    double x1, x2, x3;
    int roots_count = gsl_poly_solve_cubic(0, -b, -a, &x1, &x2, &x3);

    std::vector<EquilibriumPoint> points;
    
    std::cout << "\n--- ШАГ 1: АНАЛИЗ РАВНОВЕСИЙ ---\n";
    if (roots_count >= 1) points.push_back(analyze_point(b, x1));
    if (roots_count == 3) {
        points.push_back(analyze_point(b, x2));
        points.push_back(analyze_point(b, x3));
    }

    // Вывод информации в консоль
    for (const auto& p : points) {
        std::cout << "Точка x = " << p.value << " -> " << p.type << "\n";
    }

    // 3. Симуляция через Boost
    std::cout << "\n--- ШАГ 2: СИМУЛЯЦИЯ ДИНАМИКИ ---\n";
    BistableSystem sys(a, b);
    std::vector<double> state = { x0 };

    // -------------------------------------------------------------------------------------------------------------------------------
    // 8. Расчет состояния системы в заданные моменты времени
    // Интегрируем от t=0 до t=20 с шагом 0.1
    // integrate сама вызывает нашу систему и наблюдателя
    integrate(sys, state, 0.0, t_end, 0.1, Observer(times, x_vals));
    std::cout << "Расчет завершен. Конечное состояние: " << state[0] << "\n";

    // 4. Визуализация
    std::cout << "\n--- ШАГ 3: ПОСТРОЕНИЕ ГРАФИКА ---\n";
    
    // Настраиваем размер окна
    plt::figure_size(1000, 600);

    // Рисуем уровни равновесия (Горизонтальные линии)
    for (const auto& p : points) {
        std::vector<double> line_x(times.size(), p.value); // Линия на высоте корня
        
        // Зеленая пунктирная для устойчивых, Красная для неустойчивых
        std::string style = p.is_stable ? "g--" : "r--"; 
        
        // Используем именованные аргументы для легенды
        plt::named_plot("Равновесие " + std::to_string(p.value), times, line_x, style);
    }

    // Рисуем траекторию движения частицы (Жирная синяя линия)
    plt::named_plot("Траектория x(t)", times, x_vals, "b-");
    
    // Оформление
    plt::title("Динамика системы dx/dt = " + std::to_string(a) + " + " + std::to_string(b) + "x - x^3");
    plt::xlabel("Время (t)");
    plt::ylabel("Координата (x)");
    plt::legend(); // Показать легенду
    plt::grid(true);

    plt::save("8-time_diagram.png");
    std::cout << "График сохранен в файл 8-time_diagram.png\n";
    // -------------------------------------------------------------------------------------------------------------------------------

    // 10. Построение фазового портрета системы 
    // Подготовка данных для кривой скорости
    // Строим график для x от -2.5 до 2.5
    std::vector<double> phase_x_vals = linspace(-2.5, 2.5, 100);
    std::vector<double> phase_dxdt_vals; // Сюда запишем скорость

    for(double x : phase_x_vals) {
        // dx/dt = a + bx - x^3
        double v = a + b*x - x*x*x;
        phase_dxdt_vals.push_back(v);
    }

    // Ещё раз находим точки равновесия (через GSL), чтобы красиво их нарисовать
    int roots = gsl_poly_solve_cubic(0, -b, -a, &x1, &x2, &x3);
    
    // Векторы для хранения координат точек (отдельно устойчивые, отдельно нет)
    std::vector<double> stable_x, stable_y;
    std::vector<double> unstable_x, unstable_y;

    // Лямбда-функция для сортировки точек
    auto classify_point = [&](double root) {
        if (is_stable(b, root)) {
            stable_x.push_back(root);
            stable_y.push_back(0); // На фазовом портрете скорость равна 0
        } else {
            unstable_x.push_back(root);
            unstable_y.push_back(0);
        }
    };

    if (roots >= 1) classify_point(x1);
    if (roots == 3) {
        classify_point(x2);
        classify_point(x3);
    }

    // 3. Рисуем
    plt::figure_size(800, 600);

    // Рисуем ось X (линия y=0) черным пунктиром
    std::vector<double> zero_line(phase_x_vals.size(), 0.0);
    plt::plot(phase_x_vals, zero_line, "k--");

    // Рисуем саму функцию скорости (Синюю линию)
    plt::plot(phase_x_vals, phase_dxdt_vals, "b-");

    // Рисуем точки равновесия
    // ro = red circle (красные круги), go = green circle (зеленые круги)
    if (!unstable_x.empty()) 
        plt::plot(unstable_x, unstable_y, "ro"); 
    if (!stable_x.empty()) 
        plt::plot(stable_x, stable_y, "go");

    // Оформление
    plt::title("Фазовый портрет: Скорость dx/dt от координаты x");
    plt::xlabel("Координата x");
    plt::ylabel("Скорость dx/dt");
    plt::grid(true);

    // Добавим поясняющий текст на график
    // (Позиционирование текста примерное)
    plt::text(-2.0, 1.0, "Скорость > 0\nДвижение вправо ->");
    plt::text(1.0, -1.0, "Скорость < 0\n<- Движение влево");

    plt::save("10-phase_portrait.png");
    std::cout << "Фазовый портрет сохранен в 10-phase_portrait.png\n";
    // -------------------------------------------------------------------------------------------------------------------------------


    return 0;
}