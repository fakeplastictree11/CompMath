#include <array>
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};


unsigned int factorial(unsigned int n)
{
    if (n == 0) {return 1;}
    return n * factorial(n - 1);
}

template<typename T, unsigned int N>
std::array<T, N> Tailor(T h) {
    std::array<T, N> ans{};
    ans[0] = 1;
    for (unsigned int i = 1; i < N; i++) {
        ans[i] = pow(h, i) / factorial(i);
    }
    return ans;
}

template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
    Eigen::MatrixX<RealType> M(N + 1, N + 1);
    M.setZero();
    for (auto i = 0; i < N + 1; i++) {
        M(0, i) = 1;
    }

    for (auto j = 1; j < N + 1; j++) {
        std::array<RealType, N + 1> local_arr = Tailor<RealType, N + 1>(points[j - 1]);
        for (auto i = 1; i < N + 1; i++) {
            //for (auto k = 1; k < N + 1; k++) {std::cout << local_arr[k] << std::endl;}
            M(i, j) = local_arr[i];
        }
    }

    Eigen::VectorX<RealType> b(N + 1);
    b.setZero();
    b(L) = 1;

    Eigen::VectorX<RealType> x = M.lu().solve(b);

    DerivativeCoef<RealType, N> ans;
    ans.centralCoef = x[0];
    for (auto i = 0; i < N; i++) {
        ans.otherCoefs[i] = x[i + 1];
    }

    return ans;
}

template<typename T, unsigned int N>
std::array<T, N> create_points(T h) {
    std::array<T, N> ans;
    for (auto i = 0; i < N; i++) {ans[i] = (i + 1) * h;}
    return ans;
}

template<typename T, unsigned int N> // Создание массива со значениями функции в N точках
std::array<T, N> create_values(std::array<T, N> p) {
    std::array<T, N> v;
    for (int i = 0; i < N; i++) {
        v[i] = exp(p[i]);
    }
    return v;
}

int main() {
    std::array<double, 16> h = {1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15};
    
    /*
    // N = 3
    std::array<double, 2> p;
    std::array<double, 2> v; 
    double h_cur;
    std::array<long double, 16> err;
    for (auto i = 0; i < 16; i++) {
        h_cur = h[i];
        p = create_points<double, 2>(h_cur);
        v = create_values<double, 2>(p);

        DerivativeCoef<double, 2> coefs = calcDerivativeCoef<double, 2, 2>(p);
        
        long double ans;
        ans = coefs.centralCoef * M_E;
        for (auto j = 0; j < 2; j++) {
            ans += coefs.otherCoefs[j] * exp(1 + p[j] * h_cur);
        }

        err[i] = fabs(M_E - ans / pow(h_cur, 2));
    }
    for (auto i = 0; i < 16; i++) {std::cout << err[i] << std::endl;}

     // N = 4
    std::array<double, 3> p;
    std::array<double, 3> v; 
    double h_cur;
    std::array<long double, 16> err;
    for (auto i = 0; i < 16; i++) {
        h_cur = h[i];
        p = create_points<double, 3>(h_cur);
        v = create_values<double, 3>(p);

        DerivativeCoef<double, 3> coefs = calcDerivativeCoef<double, 3, 2>(p);
        
        long double ans;
        ans = coefs.centralCoef * M_E;
        for (auto j = 0; j < 3; j++) {
            ans += coefs.otherCoefs[j] * exp(1 + p[j] * h_cur);
        }

        err[i] = fabs(M_E - ans / pow(h_cur, 2));
    }
    for (auto i = 0; i < 16; i++) {std::cout << err[i] << std::endl;}*/

    // N = 5
    std::array<double, 4> p;
    std::array<double, 4> v; 
    double h_cur;
    std::array<long double, 16> err;
    for (auto i = 0; i < 16; i++) {
        h_cur = h[i];
        p = create_points<double, 4>(h_cur);
        v = create_values<double, 4>(p);

        DerivativeCoef<double, 4> coefs = calcDerivativeCoef<double, 4, 2>(p);
        
        long double ans;
        ans = coefs.centralCoef * M_E;
        for (auto j = 0; j < 4; j++) {
            ans += coefs.otherCoefs[j] * exp(1 + p[j] * h_cur);
        }

        err[i] = fabs(M_E - ans / pow(h_cur, 2));
    }
    for (auto i = 0; i < 16; i++) {std::cout << err[i] << std::endl;}
}
