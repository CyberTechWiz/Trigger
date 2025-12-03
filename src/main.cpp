#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
namespace odeint = boost::numeric::odeint;

//! Класс для моей системы с бистабильностью
class BistableSystem {
private:
    double a;  // параметр a
    double b;  // параметр b
    
public:
    // Конструктор
    BistableSystem(double a_val = 0.0, double b_val = 0.0) : a(a_val), b(b_val) {}
    
    // Установка параметров
    void setParameters(double a_val, double b_val) {
        a = a_val;
        b = b_val;
    }
    
    // Правая часть дифференциального уравнения: dx/dt = a + b*x - x^3
    double equation(double x) const {
        return a + b * x - x * x * x;
    }
    
    // Производная для анализа устойчивости
    double derivative(double x) const {
        return b - 3.0 * x * x;
    }
    
    // Найти корни уравнения (стационарные точки)
    std::vector<double> findEquilibriumPoints() const {
        std::vector<double> roots;
        
        // Для кубического уравнения: -x^3 + b*x + a = 0
        // Используем численное решение через GSL
        
        if (b > 0) {
            // При b > 0 могут быть 1 или 3 корня
            double discriminant = 4 * b * b * b - 27 * a * a;
            
            if (discriminant > 0) {
                // Три действительных корня
                double step = 5.0;
                double x_start = -step;
                
                for (int i = 0; i < 1000; i++) {
                    double x1 = x_start;
                    double x2 = x_start + 0.1;
                    
                    if (equation(x1) * equation(x2) < 0) {
                        // Корень на интервале [x1, x2]
                        double root = findRootBisection(x1, x2);
                        if (root != root) continue; // проверка на NaN
                        
                        // Проверяем, нет ли уже такого корня
                        bool duplicate = false;
                        for (double r : roots) {
                            if (fabs(r - root) < 1e-6) {
                                duplicate = true;
                                break;
                            }
                        }
                        if (!duplicate) {
                            roots.push_back(root);
                        }
                    }
                    x_start += 0.1;
                }
            } else {
                // Один действительный корень
                // Используем простой поиск
                double x = 0.0;
                for (int i = 0; i < 100; i++) {
                    x = x - equation(x) / derivative(x);
                }
                roots.push_back(x);
            }
        } else {
            // При b <= 0 всегда один корень
            double x = 0.0;
            for (int i = 0; i < 100; i++) {
                x = x - equation(x) / derivative(x);
            }
            roots.push_back(x);
        }
        
        // Сортируем корни
        std::sort(roots.begin(), roots.end());
        return roots;
    }
    
    // Анализ устойчивости точки равновесия
    std::string analyzeStability(double x) const {
        double deriv = derivative(x);
        if (deriv < 0) {
            return "Устойчивая";
        } else if (deriv > 0) {
            return "Неустойчивая";
        } else {
            return "Нейтральная (особая)";
        }
    }
    
    // Метод бисекции для нахождения корня
    double findRootBisection(double left, double right) const {
        double f_left = equation(left);
        double f_right = equation(right);
        
        if (f_left * f_right > 0) {
            return NAN; // Нет корня на интервале
        }
        
        for (int i = 0; i < 100; i++) {
            double mid = (left + right) / 2.0;
            double f_mid = equation(mid);
            
            if (fabs(f_mid) < 1e-12) {
                return mid;
            }
            
            if (f_left * f_mid < 0) {
                right = mid;
                f_right = f_mid;
            } else {
                left = mid;
                f_left = f_mid;
            }
        }
        
        return (left + right) / 2.0;
    }
    
    // Геттеры для параметров
    double getA() const { return a; }
    double getB() const { return b; }
};

// Класс для интегрирования системы
class SystemIntegrator {
private:
    const BistableSystem& system;
    
public:
    SystemIntegrator(const BistableSystem& sys) : system(sys) {}
    
    // Оператор для использования с Boost ODEINT
    void operator()(const std::vector<double>& x, std::vector<double>& dxdt, double t) const {
        dxdt[0] = system.getA() + system.getB() * x[0] - x[0] * x[0] * x[0];
    }
};

// Функция для построения графика
void plotVectorField(const BistableSystem& system, const std::string& filename = "vector_field.png") {
    const int n_points = 100;
    std::vector<double> x_vals(n_points);
    std::vector<double> y_vals(n_points);
    
    double x_min = -3.0;
    double x_max = 3.0;
    double dx = (x_max - x_min) / (n_points - 1);
    
    for (int i = 0; i < n_points; i++) {
        x_vals[i] = x_min + i * dx;
        y_vals[i] = system.equation(x_vals[i]);
    }
    
    plt::figure_size(800, 600);
    plt::plot(x_vals, y_vals, "b-");
    plt::axhline(0, 0, 1, {{"color", "black"}, {"linestyle", "--"}});
    plt::axvline(0, 0, 1, {{"color", "black"}, {"linestyle", "--"}});
    
    // Отметить точки равновесия
    auto equilibria = system.findEquilibriumPoints();
    for (double eq : equilibria) {
        plt::plot({eq}, {0}, "ro");
    }
    
    plt::title("Векторное поле системы: dx/dt = a + bx - x^3");
    plt::xlabel("x");
    plt::ylabel("dx/dt");
    plt::grid(true);
    plt::save(filename);
    plt::close();
}

int main() {
    std::cout << "=== Программа анализа бистабильной системы ===\n";
    std::cout << "Система: dx/dt = a + b*x - x^3\n\n";
    
    // Пример 1: Простая система
    {
        std::cout << "Пример 1: a = 0, b = 1\n";
        BistableSystem system(0, 1);
        
        auto equilibria = system.findEquilibriumPoints();
        std::cout << "Точки равновесия:\n";
        for (double eq : equilibria) {
            std::string stability = system.analyzeStability(eq);
            std::cout << "  x = " << eq << " (" << stability << ")\n";
        }
        
        plotVectorField(system, "example1.png");
        std::cout << "График сохранен в example1.png\n\n";
    }
    
    // Пример 2: Система с одной точкой равновесия
    {
        std::cout << "Пример 2: a = 1, b = 0\n";
        BistableSystem system(1, 0);
        
        auto equilibria = system.findEquilibriumPoints();
        std::cout << "Точки равновесия:\n";
        for (double eq : equilibria) {
            std::string stability = system.analyzeStability(eq);
            std::cout << "  x = " << eq << " (" << stability << ")\n";
        }
        
        plotVectorField(system, "example2.png");
        std::cout << "График сохранен в example2.png\n\n";
    }
    
    // Пример 3: Система с тремя точками равновесия
    {
        std::cout << "Пример 3: a = 0.1, b = 2\n";
        BistableSystem system(0.1, 2);
        
        auto equilibria = system.findEquilibriumPoints();
        std::cout << "Точки равновесия:\n";
        for (double eq : equilibria) {
            std::string stability = system.analyzeStability(eq);
            std::cout << "  x = " << eq << " (" << stability << ")\n";
        }
        
        // Интегрирование с помощью Boost ODEINT
        std::vector<double> x = {0.5};  // начальное условие
        odeint::runge_kutta4<std::vector<double>> stepper;
        
        std::vector<double> time_points;
        std::vector<double> solution_points;
        
        double dt = 0.01;
        double t = 0.0;
        for (int i = 0; i < 1000; i++) {
            time_points.push_back(t);
            solution_points.push_back(x[0]);
            
            SystemIntegrator integrator(system);
            stepper.do_step(integrator, x, t, dt);
            t += dt;
        }
        
        // Построение графика решения
        plt::figure_size(800, 600);
        plt::plot(time_points, solution_points, "b-");
        plt::title("Решение системы (a=0.1, b=2, x0=0.5)");
        plt::xlabel("Время t");
        plt::ylabel("x(t)");
        plt::grid(true);
        plt::save("solution.png");
        plt::close();
        
        std::cout << "График решения сохранен в solution.png\n";
        
        plotVectorField(system, "example3.png");
        std::cout << "Векторное поле сохранено в example3.png\n";
    }
    
    return 0;
}