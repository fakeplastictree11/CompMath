#include <vector>
#include <type_traits>
#include <iostream>
#include <array>
#include <cmath>


/** класс для работы с трехдиагональной матрицей **/
template<typename Type>
class ThreeDiagonalMatrix {
    private:
        std::vector<Type> a_;
        std::vector<Type> b_;
        std::vector<Type> c_;
    public:
        ThreeDiagonalMatrix(const std::vector<Type>& a, const std::vector<Type>& b, const std::vector<Type>& c): a_(a), b_(b), c_(c){}

        const Type& a(const std::size_t& i) const{return a_[i];}
        const Type& b(const std::size_t& i) const{return b_[i];}
        const Type& c(const std::size_t& i) const{return c_[i];}
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

/** Функция для решения методм  прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& A, const std::vector<cType>& d){
    std::vector<mType> p(d.size());
    std::vector<DivisType<cType, mType>> q(d.size());
    std::vector<DivisType<cType, mType>> x(d.size());

    p[0] = -A.c(0) / A.b(0);
    q[0] = d[0] / A.b(0);

    for(std::size_t i = 1; i < d.size() - 1; i++){
        p[i] = -(A.c(i) / (A.a(i - 1) * p[i - 1] + A.b(i)));
        q[i] = (d[i] - A.a(i - 1) * q[i - 1]) / (A.a(i - 1) * p[i - 1] + A.b(i));
    }

    q[d.size() - 1] = (d[d.size() - 1] - A.a(d.size() - 2) * q[d.size() - 2]) / (A.a(d.size() - 2) * p[d.size() - 2] + A.b(d.size() - 1));

    x[d.size() - 1] = q[d.size() - 1];

    for(std::size_t i = 1; i < d.size(); i++){
        x[d.size() - 1 - i] = p[d.size() - 1 - i] * x[d.size() - i] + q[d.size() - 1 - i];
    }

    return x;
}

/**
* xType - тип аргумента x.
* yType - тип значения функции y
*/
template<typename xType, typename yType>
class CubicSpline {
    /*** Какие-то поля ***/
    std::vector<xType> points;
    std::vector<yType> values;
    std::vector<xType> h;
    std::vector<yType> differences; //Разделенные разности
    std::vector<xType> coefs_a;
    std::vector<xType> coefs_b;
    std::vector<xType> coefs_c;
    std::vector<xType> coefs_d;
    
    public:
    CubicSpline( const std::vector<xType> &points_input,  // Значения x
                        const std::vector<yType>& values_input  // значения y
                        ) {
        points = points_input;
        values = values_input;
        std::vector<yType> working_arr = values;
        differences = values;
        for (auto i = 1; i < 3; i++) {
            for (auto j = i; j < values.size(); j++) {
                working_arr[j] = (differences[j] - differences[j - 1]) / (points[j] - points[j - i]);
            }
            differences = working_arr;
        } 

        std::vector<xType> lenghts(points.size() - 1);
        for (auto i = 0; i < points.size() - 1; i++) {
            lenghts[i] = points[i + 1] - points[i];
        }   
        h = lenghts;   

        std::vector<xType> as(values.size());
        for (auto i = 1; i < values.size(); i++) {
            as[i - 1] = values[i];
        } 
        coefs_a = as; 

        std::vector<xType> b(points.size() - 2);
        for (auto i = 0; i < b.size(); i++) {b[i] = 2;} 

        std::vector<xType> a(points.size() - 3);
        for (auto i = 0; i < a.size(); i++) {
            a[i] = h[i] / (h[i] + h[i + 1]);
        }

        std::vector<xType> c(points.size() - 3);
        for (auto i = 0; i < c.size(); i++) {
            c[i] = h[i + 1] / (h[i] + h[i + 1]);
        }

        std::vector<xType> d(points.size() - 2);
        for (auto i = 0; i < d.size(); i++) {
            d[i] = 6 * differences[i + 2];
        }

        coefs_c = solve(ThreeDiagonalMatrix<xType>(a, b, c), d);
        coefs_c.push_back(0);

        std::vector<xType> bs(points.size() - 1);
        bs[0] = coefs_c[0] * h[0] / 3 + (values[1] - values[0]) / (points[1] - points[0]);
        for (auto i = 1; i < points.size() - 1; i++) {
            bs[i] = coefs_c[i] * h[i] / 3 + coefs_c[i - 1] * h[i] / 6 + (values[i + 1] - values[i]) / (points[i + 1] - points[i]);
        }
        coefs_b = bs;

        std::vector<xType> ds(points.size() - 1);
        ds[0] = h[0] * coefs_c[0];
        for (auto i = 1; i < points.size() - 1; i++) {
            ds[i] = (coefs_c[i] - coefs_c[i - 1]) / h[i];
        }   
        coefs_d = ds;              
    }
                        
    yType interpolate(const xType& x) const noexcept {
        yType ans;
        int index;
        index = 0;
        for (auto i = 0; x > points[i]; i++) {
            index = i;
        }
        ans = coefs_a[index] + coefs_b[index] * (x - points[index + 1]) + coefs_c[index] * pow((x - points[index + 1]), 2) / 2 + coefs_d[index] * pow((x - points[index + 1]), 3) / 6;
        return ans;
    }

    int show() {
        for (int i = 0; i < h.size(); i++) {std::cout << h[i] << std::endl;}
        return 0;
    }
};

template<typename T> // Создание вектора с N равномерно распределенными узлами
std::vector<T> create_points(double start, double end, int N) {
    std::vector<T> p(N);
    double dx;
    dx = (end - start) / (N - 1);
    for (int i = 0; i < N - 1; i++) {
        p[i] = start + i * dx;
    }
    p[N - 1] = end;
    return p;
}

template<typename T> // Создание вектора со значениями экспоненты в N точках
std::vector<T> create_values(std::vector<T> p, int N) {
    std::vector<T> v(N);
    for (int i = 0; i < N; i++) {
        v[i] = exp(p[i]);
    }
    return v;
}

int main() {
    std::vector<double> err;
    for (auto N = 5; N < 320; N *= 2) {
        std::vector<double> p = create_points<double>(0, 10, N);
        std::vector<double> v = create_values<double>(p, N);
        std::vector<double> p1000 = create_points<double>(0, 10, 1000);

        CubicSpline<double, double> spline = CubicSpline<double, double>(p, v);
        double error, new_err;
        error = 0;
        for (int i = 0; i < 1000; i++) {
            new_err = fabs(exp(p1000[i]) - spline.interpolate(p1000[i]));
            if (error < new_err) {error = new_err;}
        }
        err.push_back(error);
    }
    for (auto i = 0; i < err.size(); i++) {std::cout << err[i] << std::endl;}
}
