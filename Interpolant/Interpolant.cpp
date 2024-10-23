#include <array>
#include <cmath>
#include <iostream>

/**
* xType - тип аргумента x.
* yType - тип значения функции y
* N - количество точек для интерполяции
*
* Рекомедую обратить внимание. Разность (xType - xType) не обязана быть типом xType
*/
template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
    /*** Какие-то поля ***/
    std::array<xType, N> points; //Узлы интерполяции
    std::array<yType, N> differences; //Разделенные разности

    public:
    //Конструктор
    NewtonInterpolator(const std::array<xType, N>& points_input, const std::array<yType, N>& values) noexcept {
        points = points_input;
        std::array<yType, N> working_arr = values;
        differences = values;
        for (auto i = 1; i < N; i++) {
            for (auto j = i; j < N; j++) {
                working_arr[j] = (differences[j] - differences[j - 1]) / (points[j] - points[j - i]);
            }
            differences = working_arr;
        }
    }

    //Интерполяция
    yType interpolate(const xType& x) const noexcept {
        yType ans, multiplication;
        ans = 0;
        for (int i = 0; i < N; i++) {
            multiplication = 1;
            for (int j = 0; j < i; j++) {
                multiplication *= x - points[j];
            }
            ans += differences[i] * multiplication;
        }
        return ans;
    }
};

template<unsigned int N> //Создание массива с N равномерно распределенными узлами
std::array<double, N> create_points(double start, double end) {
    std::array<double, N> p;
    double dx;
    dx = (end - start) / (N - 1);
    for (int i = 0; i < N - 1; i++) {
        p[i] = start + i * dx;
    }
    p[N - 1] = end;
    return p;
}

template<unsigned int N> //Создание массива с N узлами, совпадающими с корнями полинома Чебышева
std::array<double, N> create_roots(double start, double end) {
    std::array<double, N> r;
    double const1, const2;
    const1 = (end + start) / 2;
    const2 = (end - start) / 2;
    for (int i = 0; i < N; i++) {
        r[i] = const1 + const2 * cos((2 * i + 1) * M_PI / (2 * N));
    }
    return r;
}

template<unsigned int N> // Создание массива со значениями функции в N точках
std::array<double, N> create_values(std::array<double, N> p) {
    std::array<double, N> v;
    for (int i = 0; i < N; i++) {
        v[i] = exp(p[i]);
    }
    return v;
}

int main() {
    std::array<double, 1000> p1000{};
    double start = 0;
    std::array<double, 6> end{{2, 1, 0.5, 0.25, 0.125, 0.0625}};

    // Узлы распределены равномерно

    // N = 3
    std::array<double, 6> err3{};
    std::array<double, 3> p3{};
    std::array<double, 3> v3{};
    for (int i = 0; i < 6; i++){
        p3 = create_points<3>(0, end[i]);
        v3 = create_values<3>(p3);
        NewtonInterpolator<double, double, 3> interpolator3 = NewtonInterpolator<double, double, 3>(p3, v3); //Создали интерполянт

        p1000 = create_points<1000>(0, end[i]); //Считаем ошибку
        double err, new_err;
        err = 0;
        for (int i = 0; i < 1000; i++) {
            new_err = fabs(exp(p1000[i]) - interpolator3.interpolate(p1000[i]));
            if (err < new_err) {err = new_err;}
        }
        err3[i] = err;
    }

    for(int i = 0; i < 6; i++) {std::cout << err3[i] << std::endl;}
    std::cout << std::endl;

    // N = 4
    std::array<double, 6> err4{};
    std::array<double, 4> p4{};
    std::array<double, 4> v4{};
    for (int i = 0; i < 6; i++){
        p4 = create_points<4>(0, end[i]);
        v4 = create_values<4>(p4);
        NewtonInterpolator<double, double, 4> interpolator4 = NewtonInterpolator<double, double, 4>(p4, v4); //Создали интерполянт

        p1000 = create_points<1000>(0, end[i]); //Считаем ошибку
        double err, new_err;
        err = 0;
        for (int i = 0; i < 1000; i++) {
            new_err = fabs(exp(p1000[i]) - interpolator4.interpolate(p1000[i]));
            if (err < new_err) {err = new_err;}
        }
        err4[i] = err;
    }

    for(int i = 0; i < 6; i++) {std::cout << err4[i] << std::endl;}
    std::cout << std::endl;

    // N = 5
    std::array<double, 6> err5{};
    std::array<double, 5> p5{};
    std::array<double, 5> v5{};
    for (int i = 0; i < 6; i++){
        p5 = create_points<5>(0, end[i]);
        v5 = create_values<5>(p5);
        NewtonInterpolator<double, double, 5> interpolator5 = NewtonInterpolator<double, double, 5>(p5, v5); //Создали интерполянт

        p1000 = create_points<1000>(0, end[i]); //Считаем ошибку
        double err, new_err;
        err = 0;
        for (int i = 0; i < 1000; i++) {
            new_err = fabs(exp(p1000[i]) - interpolator5.interpolate(p1000[i]));
            if (err < new_err) {err = new_err;}
        }
        err5[i] = err;
    }

    for(int i = 0; i < 6; i++) {std::cout << err5[i] << std::endl;}
    std::cout << std::endl;

    // Узлы определяются как корни полинома Чебышева

    // N = 3
    std::array<double, 6> err3_ch{};
    std::array<double, 3> p3_ch{};
    std::array<double, 3> v3_ch{};
    for (int i = 0; i < 6; i++){
        p3_ch = create_roots<3>(0, end[i]);
        v3_ch = create_values<3>(p3_ch);
        NewtonInterpolator<double, double, 3> interpolator3_ch = NewtonInterpolator<double, double, 3>(p3_ch, v3_ch); //Создали интерполянт

        p1000 = create_points<1000>(0, end[i]); //Считаем ошибку
        double err, new_err;
        err = 0;
        for (int i = 0; i < 1000; i++) {
            new_err = fabs(exp(p1000[i]) - interpolator3_ch.interpolate(p1000[i]));
            if (err < new_err) {err = new_err;}
        }
        err3_ch[i] = err;
    }

    for(int i = 0; i < 6; i++) {std::cout << err3_ch[i] << std::endl;}
    std::cout << std::endl;

    // N = 4
    std::array<double, 6> err4_ch{};
    std::array<double, 4> p4_ch{};
    std::array<double, 4> v4_ch{};
    for (int i = 0; i < 6; i++){
        p4_ch = create_roots<4>(0, end[i]);
        v4_ch = create_values<4>(p4_ch);
        NewtonInterpolator<double, double, 4> interpolator4_ch = NewtonInterpolator<double, double, 4>(p4_ch, v4_ch); //Создали интерполянт

        p1000 = create_points<1000>(0, end[i]); //Считаем ошибку
        double err, new_err;
        err = 0;
        for (int i = 0; i < 1000; i++) {
            new_err = fabs(exp(p1000[i]) - interpolator4_ch.interpolate(p1000[i]));
            if (err < new_err) {err = new_err;}
        }
        err4_ch[i] = err;
    }

    for(int i = 0; i < 6; i++) {std::cout << err4_ch[i] << std::endl;}
    std::cout << std::endl;

    // N = 5
    std::array<double, 6> err5_ch{};
    std::array<double, 5> p5_ch{};
    std::array<double, 5> v5_ch{};
    for (int i = 0; i < 6; i++){
        p5_ch = create_roots<5>(0, end[i]);
        v5_ch = create_values<5>(p5_ch);
        NewtonInterpolator<double, double, 5> interpolator5_ch = NewtonInterpolator<double, double, 5>(p5_ch, v5_ch); //Создали интерполянт

        p1000 = create_points<1000>(0, end[i]); //Считаем ошибку
        double err, new_err;
        err = 0;
        for (int i = 0; i < 1000; i++) {
            new_err = fabs(exp(p1000[i]) - interpolator5_ch.interpolate(p1000[i]));
            if (err < new_err) {err = new_err;}
        }
        err5_ch[i] = err;
    }

    for(int i = 0; i < 6; i++) {std::cout << err5_ch[i] << std::endl;}
    std::cout << std::endl;
}
