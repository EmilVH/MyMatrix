//
// Created by vahna on 25.02.2020.
//

#ifndef MATRIX_ZP_FINITE_H
#define MATRIX_ZP_FINITE_H
namespace compile_time_checker{
    template<int V>
    struct get_not_zero {
        void prime() {

        }
    };

    template<>
    struct get_not_zero<0> {

    };

    template<int N, int P>
    struct check_ {
        check_<N - 1, P> x;
        get_not_zero<P - (P / N) * N> z;

        check_() {
            z.prime();
        };
    };

    template<int P>
    struct check_<1, P> {
    };

    template<int P>
    struct isPrime {
        isPrime() {
            check_<P - 1,P> x;
        }

    };
}

template<unsigned int P>
class Finite {
private:
    int module_ = P;
    int value_ = 0;

    int normalize(int value) const{
        return ((value % module_) + module_) % module_;
    }

public:
    explicit Finite<P>(int value) {
        value_ = normalize(value);
    }

    Finite<P> &operator+=(const Finite &a);

    Finite<P> &operator-=(const Finite &a);

    Finite<P> &operator*=(const Finite &a);

    Finite<P> &operator/=(const Finite &a);

    Finite<P> operator+(const Finite &a) const;

    Finite<P> operator-(const Finite &a) const;

    Finite<P> operator*(const Finite &a) const;

    Finite<P> operator/(const Finite &a) const;


    Finite<P> pow(int n) const;

    Finite<P> rev();

};

template<unsigned int P>
Finite<P> &Finite<P>::operator+=(const Finite<P> &a) {
    value_ += a.value_;
    value_ %= module_;
    return *this;
}

template<unsigned int P>
Finite<P> &Finite<P>::operator-=(const Finite &a) {
    value_ += a.value_;
    value_ = normalize(value_);
    return *this;
}

template<unsigned int P>
Finite<P> &Finite<P>::operator*=(const Finite &a) {
    value_ *= a.value_;
    value_ = normalize(value_);
    return *this;
}

template<unsigned int P>
Finite<P> &Finite<P>::operator/=(const Finite &a) {
    compile_time_checker::isPrime<P> check;
    Finite res = *this;
    res = res * a.pow(module_ - 2);
    *this = res;
    return *this;

}

template<unsigned int P>
Finite<P> Finite<P>::pow(int n) const {
    int res = 1;
    int x = value_;
    while (n > 0) {
        if (n % 2 == 1) {
            res *= x;
            --n;
            res = normalize(res);
        } else {
            x *= x;
            x = normalize(x);
            n /= 2;
        }
    }
    res = normalize(res);
    Finite ret = *this;
    ret.value_ = res;
    return ret;
}

template<unsigned int P>
Finite<P> Finite<P>::rev() {
    Finite<P> one = *this;
    one.value_ = 1;
    return one /= (*this);
}

template<unsigned int P>
Finite<P> Finite<P>::operator-(const Finite &a) const {
    Finite res = *this;
    res -= a;
    return res;
}

template<unsigned int P>
Finite<P> Finite<P>::operator+(const Finite &a) const {
    Finite res = *this;
    res += a;
    return res;
}

template<unsigned int P>
Finite<P> Finite<P>::operator*(const Finite &a) const {
    Finite res = *this;
    res *= a;
    return res;
}

template<unsigned int P>
Finite<P> Finite<P>::operator/(const Finite &a) const {
    Finite res = *this;
    res /= a;
    return res;
}


#endif //MATRIX_ZP_FINITE_H
