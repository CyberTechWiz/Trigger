#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_math.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    std::cout << "Тестирование установленных библиотек:\n";
    
    // 1. Тест Eigen
    Eigen::Matrix2d A;
    A << 1, 2, 3, 4;
    std::cout << "1. Eigen работает. Детерминант: " << A.determinant() << "\n";
    
    // 2. Тест Boost (простой пример)
    using namespace boost::numeric::odeint;
    std::vector<double> x = {1.0, 0.0};
    std::cout << "2. Boost odeint доступен\n";
    
    // 3. Тест GSL
    double result = gsl_pow_int(2, 3);
    std::cout << "3. GSL работает. 2^3 = " << result << "\n";
    
    // 4. Тест matplotlib-cpp
    std::vector<double> x_vals = {1, 2, 3, 4};
    std::vector<double> y_vals = {1, 4, 9, 16};
    
    plt::plot(x_vals, y_vals);
    plt::title("Тест matplotlib-cpp");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::save("test_plot.png");
    std::cout << "4. matplotlib-cpp работает. График сохранен в test_plot.png\n";
    
    return 0;
}