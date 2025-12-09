#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <boost/numeric/odeint.hpp>
#include "matplotlibcpp.h"
#include <gsl/gsl_poly.h> // Решатель полиномов из GSL
#include <boost/numeric/odeint.hpp> // Библиотеки Boost
#include "matplotlibcpp.h" // Библиотека рисования
#include <random> // для шума

namespace plt = matplotlibcpp;
using namespace boost::numeric::odeint;

// --- ЧАСТЬ 1: ОПИСАНИЕ СИСТЕМЫ ---
// Это "закон природы" для нашей программы.
// Компьютер будет вызывать этот код тысячи раз, чтобы узнать "куда двигаться дальше".
class BistableSystem
{
public:
    double a; // Параметр сдвига
    double b; // Параметр бифуркации
    double noise_strength; // Сила шума

    std::mt19937 gen;              // Генератор
    std::normal_distribution<> dist; // Распределение

    // Конструктор: задаем параметры один раз при создании
    BistableSystem(double param_a, double param_b, double noise) : a(param_a), b(param_b), noise_strength(noise),gen(std::random_device{}()), dist(0.0, 1.0) {}

    // Оператор (): Это стандарт Boost.
    // x    - вход: текущее значение x (где мы сейчас)
    // dxdt - выход: сюда мы записываем скорость (dx/dt)
    // t    - время (в этом уравнении оно явно не участвует, но нужно для стандарта)
    void operator()(const std::vector<double> &x, std::vector<double> &dxdt, double t)
    {
        // Наше уравнение: dx/dt = a + bx - x^3
        // x[0], потому что вектор может хранить много переменных, но у нас одна.
        double deterministic = a + b * x[0] - (x[0] * x[0] * x[0]);
        double random_kick = noise_strength * dist(gen);
        dxdt[0] = deterministic + random_kick;
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
    double noise_val = 0.5;      // Сила шума (Пункт 12)
    BistableSystem system(a, b, noise_val);

// 3. Задание начальных условий
    // Начнем с x = 0.1. Система должна "скатиться" в устойчивое состояние
    std::vector<double> x_state = {0.1};

    // Векторы, куда будем писать историю
    std::vector<double> times;
    std::vector<double> x_vals;

// -------------------------------------------------------------------------------------------------------------------------------
// 2. Определение точек равновесия и их типа
    // --- ПОДГОТОВКА ДЛЯ GSL ---
    // GSL решает x^3 + c1*x^2 + c2*x + c3 = 0
    // Наше уравнение: x^3 - b*x - a = 0
    double x1, x2, x3;
    int roots_count = gsl_poly_solve_cubic(0, -b, -a, &x1, &x2, &x3);

    std::vector<EquilibriumPoint> points;
    
    std::cout << "\n1. Определение положений равновесий и их типов\n";
    if (roots_count >= 1) points.push_back(analyze_point(b, x1));
    if (roots_count == 3) {
        points.push_back(analyze_point(b, x2));
        points.push_back(analyze_point(b, x3));
    }

    // Вывод информации в консоль
    for (const auto& p : points) {
        std::cout << "Точка x = " << p.value << " -> " << p.type << "\n";
    }

    // Симуляция через Boost
    std::vector<double> state = { x0 };
// -------------------------------------------------------------------------------------------------------------------------------

// 4, 5. Задание возмущающих и управляющих воздействий
    std::cout << "\n4,5. Задание возмущабщих и управляющих воздействий\n";
    // Параметры управления (изменение параметра a)
    double t_control_start = 5.0; 
    double t_control_end = 10.0;
    double a_modified = 2.0;     // Значение 'a' во время управления (наклон ямы)

    // Параметры возмущения (резкий удар по x)
    double t_perturbation = 15.0;
    double kick_size = -2.0;     // Насколько сильно пнуть шарик
    bool is_kicked = false;      // Флаг, чтобы пнуть только один раз
// -------------------------------------------------------------------------------------------------------------------------------

// 6. РАСЧЕТ ХАРАКТЕРИСТИК ВЕКТОРНОГО ПОЛЯ СКОРОСТЕЙ
    std::cout << "\n6. Расчет характеристик векторного поля скоростей\n";
    // Анализируем поле на интервале x от -2.0 до 2.0
    double x_field_start = -2.0;
    double x_field_end = 2.0;
    int field_steps = 10; // Сколько точек проверить
    double step_val = (x_field_end - x_field_start) / field_steps;
    
    std::cout << "--------------------------------------------------\n";
    std::cout << std::setw(10) << "X" << " | " << std::setw(15) << "Скорость (dx/dt)" << " | " << "Направление" << "\n";
    std::cout << "--------------------------------------------------\n";

    double max_speed = 0.0;

    for(int i = 0; i <= field_steps; ++i) {
        double current_x = x_field_start + i * step_val;
        // Считаем скорость в этой точке: v = a + bx - x^3
        // Используем параметры a и b, заданные в начале
        double v = a + b * current_x - pow(current_x, 3);
        
        // Определяем направление
        std::string direction = (v > 0) ? "-> (Вправо)" : (v < 0) ? "<- (Влево)" : "STOP";
        if(std::abs(v) < 1e-6) direction = "Равновесие";

        // Поиск максимальной скорости (по модулю) для статистики
        if(std::abs(v) > max_speed) max_speed = std::abs(v);

        std::cout << std::fixed << std::setprecision(2) << std::setw(10) << current_x 
                  << " | " << std::setw(15) << v 
                  << " | " << direction << "\n";
    }
    std::cout << "Максимальная скорость в диапазоне: " << max_speed << "\n";
// -------------------------------------------------------------------------------------------------------------------------------

// 7, 12. МОДЕЛИРОВАНИЕ С ШУМОМ, УПРАВЛЕНИЕМ И ВОЗМУЩЕНИЯМИ
    std::cout << "\n7. Моделирование управляющих воздействий и возмущающих воздействий (шум)...\n";
    
    // Создаем систему с шумом
    BistableSystem sys(a, b, noise_val);
    
    std::vector<double> x = { x0 }; // Текущее состояние
    
    // Очищаем векторы для новых данных
    times.clear(); 
    x_vals.clear();
    std::vector<double> a_vals; // Для графика управления

    // Используем ручной цикл по времени вместо integrate(), 
    // чтобы менять параметры внутри цикла
    double dt = 0.05;
    runge_kutta4<std::vector<double>> stepper; // Шаговик

    for (double t = 0; t <= t_end; t += dt) {
        
        // --- Реализация Управляющего воздействия ---
        // Если время внутри интервала управления, меняем параметр системы
        if (t >= t_control_start && t < t_control_end) {
            sys.a = a_modified; 
        } else {
            sys.a = a; // Возвращаем исходное значение
        }

        // --- Реализация Возмущающего воздействия ---
        // Если пришло время удара
        if (t >= t_perturbation && !is_kicked) {
            x[0] += kick_size; // Мгновенно меняем координату
            is_kicked = true;
            std::cout << "--> УДАР (Возмущение) в момент t=" << t << "\n";
        }

        // Запись данных
        times.push_back(t);
        x_vals.push_back(x[0]);
        a_vals.push_back(sys.a); // Записываем, какой был параметр a

        // Сделать шаг физики
        stepper.do_step(sys, x, t, dt);
    }

    // --- ВИЗУАЛИЗАЦИЯ СЛОЖНОЙ ДИНАМИКИ ---
    // --- ИСПРАВЛЕННАЯ ВИЗУАЛИЗАЦИЯ СЛОЖНОЙ ДИНАМИКИ ---
    // Вместо subplot мы сохраним два красивых отдельных файла.
    
    // Рисуем график состояния X(t)
    plt::clf(); // <--- ВАЖНО: Очищаем память от предыдущего графика (Бифуркации)
    plt::figure_size(1000, 600);
    
    plt::plot(times, x_vals, "b-");
    plt::title("Dynamics with Control and Noise");
    plt::ylabel("Coordinate X");
    plt::xlabel("Time t");
    plt::grid(true);
    
    // Добавим подписи событий, если хотите
    plt::text(t_control_start, -0.5, "Управление");
    plt::text(t_perturbation, x_vals.back(), "Удар"); // Примерная позиция

    plt::save("7-dynamics.png");
    std::cout << "График динамики сохранен в 7-dynamics.png\n";

    // 2. Рисуем график управления A(t)
    plt::clf(); // <--- ВАЖНО: Снова очищаем холст для нового рисунка
    plt::figure_size(1000, 400); // Можно сделать его поменьше
    
    plt::plot(times, a_vals, "r-");
    plt::title("Control Parameter a(t)");
    plt::ylabel("Parameter a");
    plt::xlabel("Time t");
    plt::grid(true);

    plt::save("7-control.png");
    std::cout << "График управления сохранен в 7-control.png\n";
// -------------------------------------------------------------------------------------------------------------------------------

// 8. Расчет состояния системы в заданные моменты времени
    // Интегрируем от t=0 до t=20 с шагом 0.1
    // integrate сама вызывает нашу систему и наблюдателя
    std::cout << "\n8. Расчет состояния системы в заданные моменты времени \n";
    integrate(system, state, 0.0, t_end, 0.1, Observer(times, x_vals));
    std::cout << "Расчет завершен. Конечное состояние: " << state[0] << "\n";

    // Визуализация
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

    // 9. ВЫВОД РАССЧИТАННЫХ ЗНАЧЕНИЙ СОСТОЯНИЯ В ТАБЛИЧНОМ ВИДЕ
    std::cout << "\n9. Табличный вывод результатов симуляции (первые 10 шагов и последние 5)...\n";
    std::cout << "-----------------------------------------\n";
    std::cout << std::setw(10) << "Время (t)" << " | " << std::setw(15) << "Состояние (x)" << " | " << "Скорость (dx/dt)\n";
    std::cout << "-----------------------------------------\n";

    // Поскольку точек очень много (20 / 0.1 = 200 точек), выведем не всё, чтобы не засорять консоль
    size_t total_points = times.size();
    
    for (size_t i = 0; i < total_points; ++i) {
        // Выводим первые 10 точек ИЛИ последние 5 точек
        if (i < 10 || i > total_points - 6) {
            double t_curr = times[i];
            double x_curr = x_vals[i];
            // Рассчитаем мгновенную скорость для справки
            double v_curr = a + b * x_curr - pow(x_curr, 3);

            std::cout << std::fixed << std::setprecision(4) << std::setw(10) << t_curr 
                      << " | " << std::setw(15) << x_curr 
                      << " | " << v_curr << "\n";
        }
        // Если это 10-я точка, поставим многоточие
        if (i == 10) {
            std::cout << "   ...    (данные скрыты)    ...   \n";
        }
    }
    std::cout << "-----------------------------------------\n";

// -------------------------------------------------------------------------------------------------------------------------------

// 10. Построение фазового портрета системы 
    // Подготовка данных для кривой скорости
    // Строим график для x от -2.5 до 2.5
    std::cout << "\n10. Построение фазового портрета\n";
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

    // Рисуем
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
    // plt::text(-2.0, 1.0, "Скорость > 0\nДвижение вправо ->");
    // plt::text(1.0, -1.0, "Скорость < 0\n<- Движение влево");

    plt::save("10-phase_portrait.png");
    std::cout << "Фазовый портрет сохранен в 10-phase_portrait.png\n";
// -------------------------------------------------------------------------------------------------------------------------------

// 11. ПОСТРОЕНИЕ ПАРАМЕТРИЧЕСКОГО ПОРТРЕТА (БИФУРКАЦИОННАЯ ДИАГРАММА)
std::cout << "\n11. Построение параметрического портрета...\n";
    
    // ИСПРАВЛЕНИЕ: Создаем ПАРЫ векторов (X, Y) для каждого типа точек
    std::vector<double> b_stable, x_stable;     // Для зеленых точек
    std::vector<double> b_unstable, x_unstable; // Для красных точек
    
    // Пробегаем параметром b от -2 до 5
    for (double b_param = -2.0; b_param <= 5.0; b_param += 0.05) {
        double r1, r2, r3;
        // Решаем уравнение для текущего b (при фиксированном начальном a)
        int n = gsl_poly_solve_cubic(0, -b_param, -a, &r1, &r2, &r3);
        
        // Лямбда-функция теперь добавляет данные в СООТВЕТСТВУЮЩИЕ пары векторов
        auto save_if_root = [&](double root) {
            // Используем вашу функцию analyze_point или is_stable
            // Здесь для простоты используем проверку лямбды напрямую или вызовем analyze_point
            EquilibriumPoint pt = analyze_point(b_param, root); // Используем вашу функцию
            
            if (pt.is_stable) {
                b_stable.push_back(b_param);   // X координата для устойчивых
                x_stable.push_back(root);      // Y координата для устойчивых
            } else {
                b_unstable.push_back(b_param); // X координата для неустойчивых
                x_unstable.push_back(root);    // Y координата для неустойчивых
            }
        };

        if (n >= 1) save_if_root(r1);
        if (n == 3) { save_if_root(r2); save_if_root(r3); }
    }

    plt::figure_size(800, 600);
    
    // Теперь передаем согласованные векторы: (b_stable, x_stable)
    if (!b_stable.empty())
        plt::scatter(b_stable, x_stable, 10.0, {{"color", "green"}, {"label", "Устойчивые"}});
    
    // И (b_unstable, x_unstable)
    if (!b_unstable.empty())
        plt::scatter(b_unstable, x_unstable, 10.0, {{"color", "red"}, {"label", "Неустойчивые/Полууст."}});
        
    plt::title("11. Параметрический портрет (Бифуркационная диаграмма)");
    plt::xlabel("Параметр b");
    plt::ylabel("Положения равновесия x");
    plt::grid(true);
    // plt::legend(); // Легенда в scatter иногда работает некорректно в старых версиях, если будет ошибка - закомментируйте
    plt::save("11-bifurcation.png");
    std::cout << "Параметрический портрет сохранен в 11-bifurcation.png\n";
// -------------------------------------------------------------------------------------------------------------------------------



    return 0;
}