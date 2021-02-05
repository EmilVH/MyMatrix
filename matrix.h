//
// Created by vahna on 09.03.2020.
//

#ifndef MATRIX_ZP_MATRIX_H
#define MATRIX_ZP_MATRIX_H

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

namespace compile_time_checker {
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
        get_not_zero<P - (P / N) * N> z;
        check_<(N + 1) * (N * N < P), P> x;

        check_() {
            z.prime();
        };
    };

    template<int P>
    struct check_<1, P> {
    };

    template<int P>
    struct check_<0, P> {
    };

    template<int P>
    struct isPrime {
        isPrime() {
            check_<2, P> x;
        }

    };
}

//Begin of fast rational
/*
typedef std::complex<long double> cld;

const long double PI = acosl(-1.0);

size_t findUpperDegreeOfTwo(size_t v) {
    size_t n = 1;
    while (n < v)
        n *= 2;
    return n;
}

size_t reverseBits(size_t i, size_t n) {
    size_t ans = 0;
    for (; n > 1; n /= 2) {
        ans = ans * 2 + (i % 2);
        i /= 2;
    }
    return ans;
}

typedef std::vector<cld>::iterator TIt;

void makeGeneralFFT(TIt begin, TIt end, cld q) {
    size_t n = end - begin;
    if (n == 1)
        return;

    for (size_t i = 0; i < n; ++i) {
        size_t revi = reverseBits(i, n);
        if (i < revi)
            std::swap(*(begin + i), *(begin + revi));
    }

    for (size_t l = 2; l <= n; l *= 2) {
        cld ql = q;
        for (size_t ll = n; ll > l; ll /= 2)
            ql *= ql;


        for (auto itbegin = begin; itbegin != end; itbegin += l) {
            cld w(1.0, 0.0);
            auto lit = itbegin;
            auto rit = itbegin + l / 2;
            auto eit = itbegin + l;
            for (; rit < eit; ++lit, ++rit) {
                cld u = *lit;
                cld v = w * (*rit);
                *lit = u + v;
                *rit = u - v;
                w *= ql;
            }
        }
    }
}

std::vector<cld> makeGeneralFFT(std::vector<cld> a, cld q) {
    makeGeneralFFT(a.begin(), a.end(), q);
    return a;
}


std::vector<cld> makeFFT(const std::vector<cld> &a) {
    long double ang = 2.0 * PI / a.size();
    return makeGeneralFFT(a, cld(cosl(ang), sinl(ang)));
}

std::vector<cld> makeInverseFFT(std::vector<cld> a) {
    long double ang = 2.0 * PI / a.size();
    a = makeGeneralFFT(a, cld(cosl(ang), -sinl(ang)));
    for (size_t i = 0; i < a.size(); ++i)
        a[i] /= a.size();
    return a;
}


cld makeComplexNumber(int x) {
    return cld(x, 0.0);
}


std::vector<cld> makeComplexVector(const std::vector<int> &a, size_t n) {
    std::vector<cld> ca(n, cld(0.0, 0.0));
    for (size_t i = 0; i < a.size(); ++i)
        ca[i] = makeComplexNumber(a[i]);
    return ca;
}


std::vector<int> makeTVector(const std::vector<cld> &a) {
    std::vector<int> ans(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        ans[i] = static_cast<int>(floorl(a[i].real() + 0.5));
    return ans;
}

template<typename T>
std::vector<T> multiplicatePolynoms(const std::vector<T> &a,
                                    const std::vector<T> &b) {
    size_t n = findUpperDegreeOfTwo(std::max(a.size(), b.size()));
    n *= 2;

    std::vector<cld> ca = makeComplexVector(a, n);
    std::vector<cld> cb = makeComplexVector(b, n);

    ca = makeFFT(ca);
    cb = makeFFT(cb);

    for (size_t i = 0; i < n; ++i)
        ca[i] *= cb[i];

    std::vector<cld> cc = makeInverseFFT(ca);
    return makeTVector(cc);
}

template<class curr>
void my_reverse(curr first, curr last) {
    while ((first != last) && (first != --last)) {
        std::swap(*first++, *last);
    }
}

struct BigInteger {
public:


    BigInteger();

    BigInteger(int a);

    BigInteger(const std::string &s);

    explicit operator bool() const;

    explicit operator int() const;

    std::string toString() const;

    int compare(const BigInteger &a, const BigInteger &b, bool byAbs) const;

    BigInteger abs() const;

    BigInteger &operator+=(const BigInteger &a);

    BigInteger &operator-=(const BigInteger &a);

    BigInteger &operator*=(const BigInteger &a);

    BigInteger &operator/=(const BigInteger &b);

    BigInteger &operator%=(const BigInteger &b);

    friend bool operator>(const BigInteger &a, const BigInteger &b);

    friend bool operator<(const BigInteger &a, const BigInteger &b);

    friend bool operator==(const BigInteger &a, const BigInteger &b);

    friend bool operator!=(const BigInteger &a, const BigInteger &b);

    friend bool operator<=(const BigInteger &a, const BigInteger &b);

    friend bool operator>=(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator+(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator-(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator*(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator/(const BigInteger &a, const BigInteger &b);

    friend BigInteger operator%(const BigInteger &a, const BigInteger &b);

    BigInteger &operator++();

    BigInteger &operator--();

    BigInteger operator++(int);

    BigInteger operator--(int);

    BigInteger operator-() const;


public:
    std::vector<int> numb_;

    bool sign_ = true;

    static const int BASE_ = static_cast<int>(10);

    void deleteLeadingZeroes_();

    void addAbs_(BigInteger &a, const BigInteger &b);

    void substractAbs_(const BigInteger &first, const BigInteger &second);

    void div_(const BigInteger &a, const BigInteger &b, BigInteger &res, BigInteger &mod);

    int getRes_(int a, int b) const;

    void normalize_();
};

BigInteger::BigInteger() {
    numb_.push_back(0);
}

BigInteger::BigInteger(int a) {
    if (a < 0) {
        sign_ = false;
    }
    a = std::abs(a);
    if (a == 0) {
        numb_.push_back(0);
    } else {
        while (a > 0) {
            numb_.push_back(a % 10);
            a /= 10;
        }
    }
}

BigInteger::BigInteger(const std::string &s) {
    int strt = 0;
    if (s[0] == '-') {
        sign_ = false;
        strt++;
    }
    if (s == "0") {
        numb_.push_back(0);
    } else {
        for (int i = static_cast<int>(s.size() - 1); i >= strt; --i) {
            numb_.push_back(s[i] - '0');
        }
    }
}

void BigInteger::deleteLeadingZeroes_() {
    for (int i = static_cast<int>(numb_.size()) - 1; i >= 0; --i) {
        if (numb_[i] == 0) {
            numb_.pop_back();
        } else {
            break;
        }
    }
    if (numb_.empty()) {
        numb_.push_back(0);
    }
}

void BigInteger::addAbs_(BigInteger &a, const BigInteger &b) {
    int over = 0;
    size_t max_size = std::max(a.numb_.size(), b.numb_.size());
    for (size_t i = 0; i < max_size; ++i) {
        int first = 0;
        int second = 0;
        if (a.numb_.size() > i) first = a.numb_[i];
        if (b.numb_.size() > i) second = b.numb_[i];
        if (a.numb_.size() > i) {
            a.numb_[i] = first + second + over;
        } else {
            a.numb_.push_back(first + second + over);
        }
        over = 0;
        if (a.numb_[i] >= BASE_) {
            a.numb_[i] -= BASE_;
            over = 1;
        }
    }
    if (over > 0) {
        a.numb_.push_back(1);
    }

}

void BigInteger::substractAbs_(const BigInteger &first, const BigInteger &second) {
    int owed = 0;
    numb_.resize(std::max(first.numb_.size(), second.numb_.size()));
    for (size_t i = 0; i < second.numb_.size(); ++i) {
        if (first.numb_[i] < (second.numb_[i] + owed)) {
            numb_[i] = first.numb_[i] + BASE_ - (second.numb_[i] + owed);
            owed = 1;
        } else {
            numb_[i] = first.numb_[i] - (second.numb_[i] + owed);
            owed = 0;
        }
    }
    for (size_t i = second.numb_.size(); i < first.numb_.size(); ++i) {
        if (first.numb_[i] < owed) {
            numb_[i] = first.numb_[i] + BASE_ - (owed);
            owed = 1;
        } else {
            numb_[i] = first.numb_[i] - (owed);
            owed = 0;
        }
    }
    deleteLeadingZeroes_();

}

void BigInteger::div_(const BigInteger &a, const BigInteger &b, BigInteger &res, BigInteger &mod) {
    for (int i = static_cast<int>(a.numb_.size() - 1); i > -1; --i) {
        mod *= static_cast<int>(BASE_);
        mod += a.numb_[i];
        unsigned int l = 0;
        unsigned int r = BASE_;
        while (r - l > 1) {
            unsigned int m = (r + l) / 2;
            BigInteger tmp = b.abs() * static_cast<int>(m);
            if (tmp > mod) {
                r = m;
            } else {
                l = m;
            }
        }
        res.numb_.push_back(l);
        mod -= b.abs() * static_cast<int>(l);
    }
    my_reverse(res.numb_.begin(), res.numb_.end());
    res.sign_ = !(b.sign_ ^ a.sign_);
    mod.sign_ = !(b.sign_ ^ a.sign_);
    res.deleteLeadingZeroes_();
    mod.deleteLeadingZeroes_();
}

BigInteger::operator bool() const {
    return (*this) != 0;
}

BigInteger::operator int() const {
    int ans = 0;
    for (int i = static_cast<int>(numb_.size() - 1); i >= 0; --i) {
        ans *= BASE_;
        ans += numb_[i];
    }
    ans *= sign_;
    return ans;
}

std::string BigInteger::toString() const {
    std::string ans;
    if (!sign_ && compare(*this, 0, true) != 0) {
        ans.push_back('-');
    }
    for (size_t i = 0; i < numb_.size(); ++i) {
        ans.push_back(numb_[numb_.size() - 1 - static_cast<int>(i)] + '0');
    }
    return ans;
}

int BigInteger::compare(const BigInteger &a, const BigInteger &b, bool byAbs) const {
    if (byAbs) {
        if (a.numb_.size() != b.numb_.size()) {
            return getRes_(a.numb_.size(), b.numb_.size());
        }
        for (int i = static_cast<int>(a.numb_.size()) - 1; i >= 0; --i) {
            if (a.numb_[i] != b.numb_[i]) {
                return getRes_(a.numb_[i], b.numb_[i]);
            }
        }
        return 0;
    } else {
        if (a.sign_ == b.sign_) {
            return a.sign_ ? compare(a, b, true) : -compare(a, b, true);
        } else {
            return a.sign_ ? 1 : -1;
        }
    }
}

BigInteger BigInteger::abs() const {
    BigInteger res = *this;
    res.sign_ = true;
    return res;
}

BigInteger &BigInteger::operator+=(const BigInteger &a) {
    if (sign_ == a.sign_) {
        addAbs_(*this, a);
    } else {

        if (compare(*this, a, true) > -1) {
            substractAbs_(*this, a);
        } else {
            substractAbs_(a, *this);
            sign_ = !sign_;
        }
    }
    normalize_();
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &a) {

    if (sign_ == a.sign_) {
        if (compare(*this, a, true) > -1) {
            substractAbs_(*this, a);
        } else {
            substractAbs_(a, *this);
            sign_ = !sign_;
        }
    } else {
        addAbs_(*this, a);
    }
    normalize_();

    return *this;
}

BigInteger &BigInteger::operator*=(const BigInteger &a) {
    BigInteger res;
    res.numb_.resize(numb_.size() + a.numb_.size() + 2, 0);
    res.numb_ = multiplicatePolynoms(numb_, a.numb_);
    int k = 0;
    for (size_t i = 0; i < res.numb_.size(); i++) {
        res.numb_[i] += k;
        k = res.numb_[i] / 10;
        res.numb_[i] %= 10;
    }
    res.deleteLeadingZeroes_();
    res.sign_ = !(sign_ ^ a.sign_);
    *this = res;
    normalize_();
    return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &b) {
    BigInteger res;
    BigInteger mod;
    div_(*this, b, res, mod);
    *this = res;
    normalize_();
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &b) {
    BigInteger res;
    BigInteger mod;
    div_(*this, b, res, mod);
    *this = mod;
    normalize_();
    return *this;
}

bool operator>(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) == 1;
}

bool operator<(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) == -1;
}

bool operator==(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) == 0;
}

bool operator!=(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) != 0;
}

bool operator<=(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) < 1;
}

bool operator>=(const BigInteger &a, const BigInteger &b) {
    return a.compare(a, b, false) > -1;
}

BigInteger operator+(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c += b;
    return c;
}

BigInteger operator-(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c -= b;
    return c;
}

BigInteger operator*(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c *= b;
    return c;
}

BigInteger operator/(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c /= b;
    return c;
}

BigInteger operator%(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c %= b;
    return c;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &b) {
    return out << b.toString();
}

BigInteger &BigInteger::operator++() {
    return (*this) += 1;
}

BigInteger &BigInteger::operator--() {
    return (*this) -= 1;
}

BigInteger BigInteger::operator++(int) {
    BigInteger ans = (*this);
    ++(*this);
    return ans;
}

BigInteger BigInteger::operator-() const {
    BigInteger ans = *this;
    ans.sign_ ^= 1;
    return ans;
}

BigInteger BigInteger::operator--(int) {
    BigInteger ans = (*this);
    --(*this);
    return ans;
}

void BigInteger::normalize_() {
    if (compare(*this, 0, true) == 0) {
        sign_ = true;
    }
}

int BigInteger::getRes_(int a, int b) const {
    if (a < b) {
        return -1;
    }
    if (a > b) {
        return 1;
    }
    return 0;
}

std::istream &operator>>(std::istream &in, BigInteger &b) {
    std::string s;
    in >> s;

    b = s;
    return in;
}

BigInteger get_gcd(BigInteger a, BigInteger b) {
    while (b) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}

struct Rational {
public:

    Rational();

    Rational(long long a);

    Rational(int a);

    Rational(const BigInteger &a);

    std::string toString() const;

    explicit operator double();

    std::string asDecimal(size_t precision = 0) const;

    int compare(const Rational &a, bool byAbs) {
        return x_.compare(x_ * a.y_, a.x_ * y_, byAbs);
    }

    Rational &operator+=(const Rational &a);

    Rational &operator-=(const Rational &a);

    Rational &operator*=(const Rational &a);

    Rational &operator/=(const Rational &a);

    bool operator>(const Rational &a);

    bool operator<(const Rational &a);

    bool operator==(const Rational &a);

    bool operator!=(const Rational &a);

    bool operator<=(const Rational &a);

    bool operator>=(const Rational &a);

    friend Rational operator+(const Rational &a, const Rational &b);

    friend Rational operator-(const Rational &a, const Rational &b);

    friend Rational operator*(const Rational &a, const Rational &b);

    friend Rational operator/(const Rational &a, const Rational &b);

    friend std::istream &operator>>(std::istream &in, Rational &b) {
        std::string s;
        in >> s;

        int strt = 0;
        if (s[0] == '-') {
            b.x_.sign_ = false;
            strt++;
        }
        b.x_.numb_.resize(0);
        bool fl = false;
        if (s == "0") {
            b.x_.numb_.push_back(0);
        } else {
            for (int i = static_cast<int>(s.size() - 1); i >= strt; --i) {
                if (s[i] != '.') {
                    b.x_.numb_.push_back(s[i] - '0');
                } else {
                    fl = true;
                }
            }
        }
        b.y_ = 1;
        if (fl) {
            for (int i = static_cast<int>(s.length() - 1); i >= 0; i--) {
                if (s[i] != '.') {
                    b.y_ *= 10;
                } else {
                    break;
                }
            }
        }

        return in;
    }

    Rational operator-() const;

    BigInteger cut_value = 100000;

    void cut() {
        size_t length = 5;
        while (((std::min(x_.numb_.size(), y_.numb_.size()) > 30) ||
                std::max(x_.numb_.size(), y_.numb_.size()) > 120) &&
               (x_.numb_.size() > length && y_.numb_.size() > length)) {
            for (size_t i = 0; i + length < x_.numb_.size(); i++) {
                x_.numb_[i] = x_.numb_[i + length];
            }
            x_.numb_.resize(x_.numb_.size() - length);
            for (size_t i = 0; i + length < y_.numb_.size(); i++) {
                y_.numb_[i] = y_.numb_[i + length];
            }
            y_.numb_.resize(y_.numb_.size() - length);
            //std::cout << "cutted" << std::endl;
        }
        if (y_.numb_.size() > 200) {

            x_ = 0;
            y_ = 1;

        }
    }

private:
    BigInteger x_;
    BigInteger y_;

    void simplify_();

};

void Rational::simplify_() {
    BigInteger gcd = get_gcd(x_.abs(), y_.abs());
    x_ /= gcd;
    y_ /= gcd;

}

Rational::Rational() {
    x_ = 0;
    y_ = 1;
}

Rational::Rational(long long a) {
    x_ = a;
    y_ = 1;
}

Rational::Rational(int a) {
    x_ = a;
    y_ = 1;
}

Rational::Rational(const BigInteger &a) {
    x_ = a;
    y_ = 1;
}

std::string Rational::toString() const {

    std::string ans;
    ans += x_.toString();
    if (y_ != 1 && x_ != 0) {
        ans += "/";
        ans += y_.toString();
    }
    return ans;
}


Rational::operator double() {
    std::string s = asDecimal(50);
    double ans = std::stod(s);
    return ans;
}

Rational &Rational::operator+=(const Rational &a) {
    x_ = (x_ * a.y_) + (a.x_ * y_);
    y_ = y_ * a.y_;
    //simplify_();
    cut();
    return *this;
}

Rational &Rational::operator-=(const Rational &a) {
    x_ = (x_ * a.y_) - (a.x_ * y_);
    y_ = y_ * a.y_;
    //simplify_();
    cut();
    return *this;
}

Rational operator+(const Rational &a, const Rational &b) {
    Rational c = a;
    c += b;
    return c;
}

Rational operator-(const Rational &a, const Rational &b) {
    Rational c = a;
    c -= b;
    return c;
}

Rational operator*(const Rational &a, const Rational &b) {
    Rational c = a;
    c *= b;
    return c;
}

Rational operator/(const Rational &a, const Rational &b) {
    Rational c = a;
    c /= b;
    return c;
}

Rational Rational::operator-() const {
    Rational ans = *this;
    ans.x_ = -ans.x_;
    return ans;
}

Rational &Rational::operator/=(const Rational &a) {
    //if (&a == this) return (*this) = 1;
    x_ *= a.y_;
    y_ *= a.x_;
    if (y_ < 0) {
        y_ = y_.abs();
        x_ = -x_;
    }
    //simplify_();
    cut();
    return *this;
}

Rational &Rational::operator*=(const Rational &a) {
    x_ *= a.x_;
    y_ *= a.y_;
    //simplify_();
    cut();
    //std::cout << x_.numb_.size() << std::endl << y_.numb_.size() << std::endl;
    return *this;
}

bool Rational::operator>(const Rational &a) {
    return (*this).compare(a, false) == 1;
}

bool Rational::operator<(const Rational &a) {
    return (*this).compare(a, false) == -1;
}

bool Rational::operator==(const Rational &a) {
    Rational EPS;
    EPS.x_ = 1;
    EPS.y_ = 100000;
    return ((*this) - a).compare(EPS, true) == -1;
}

bool Rational::operator!=(const Rational &a) {
    Rational EPS;
    EPS.x_ = 1;
    EPS.y_ = BigInteger("1000000");
    return ((*this) - a).compare(EPS, true) == 1;
}

bool Rational::operator<=(const Rational &a) {
    return (*this).compare(a, false) != 1;
}

bool Rational::operator>=(const Rational &a) {
    return (*this).compare(a, false) != -1;
}

std::string Rational::asDecimal(size_t precision) const {
    BigInteger mul = 1;
    for (size_t i = 0; i < precision; i += 10) {
        mul *= BigInteger("10000000000");
    }
    //std::cout << "precision calc" << std::endl;
    BigInteger res = x_ * BigInteger(mul) / y_;
    //std::cout << "multiplicated" << std::endl;
    std::string sign_string;
    if (res < 0) {
        res = -res;
        sign_string = "-";
    }
    std::string pre_ans = (res).toString();
    my_reverse(pre_ans.begin(), pre_ans.end());
    pre_ans.resize(std::max(pre_ans.size(), precision + 1), '0');
    my_reverse(pre_ans.begin(), pre_ans.end());
    std::string ans = sign_string + pre_ans.substr(0, pre_ans.size() - precision);
    if (precision != 0) {
        return ans + "." + pre_ans.substr(pre_ans.size() - precision, pre_ans.length());
    } else {
        return ans;
    }
}
*/
//End of fast rational
template<unsigned int P>
class Finite {
private:
    unsigned int module_ = P;
public:
    int value_ = 0;

    int normalize(long long value) const {
        if (value < 0) {
            if (-value < static_cast<long long>(P)) {
                value += static_cast<long long>(P);
            } else {
                value %= static_cast<long long>(P);
                value += static_cast<long long>(P);
            }
        } else {
            if (value < static_cast<long long>(P)) {

            } else if (value < static_cast<long long>(2 * P)) {
                value -= static_cast<long long>(P);
            } else {
                value %= static_cast<long long>(P);
            }
        }
        return value;
    }

public:
    Finite<P>() {
        value_ = normalize(0);
    }

    Finite<P>(int value) {
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

    Finite<P> &operator++() {
        ++value_;
        value_ = normalize(value_);
        return (*this);
    }

    Finite<P> &operator--() {
        --value_;
        value_ = normalize(value_);
        return (*this);
    }

    bool operator>(const Finite<P> &x) {

        return (value_ > x.value_);
    }

    bool operator>=(const Finite<P> &x) {

        return (value_ >= x.value_);
    }

    bool operator<=(const Finite<P> &x) {

        return (value_ <= x.value_);
    }

    Finite<P> operator++(int) {
        Finite<P> ans = (*this);
        ++value_;
        value_ = normalize(value_);
        return ans;
    }

    Finite<P> operator--(int) {
        Finite<P> ans = (*this);
        --value_;
        value_ = normalize(value_);
        return ans;
    }

    Finite<P> operator-() {
        Finite<P> ans = *this;
        ans.value_ = normalize(ans.value_ * -1);
        return ans;
    }

    bool operator==(const int &other) const {
        return normalize(value_) == normalize(other);
    }

    bool operator!=(const int &other) const {
        return normalize(value_) != normalize(other);
    }

    bool operator==(const Finite<P> &other) const {
        return normalize(value_) == normalize(other.value_);
    }

    bool operator!=(const Finite<P> &other) const {

        return normalize(value_) != normalize(other.value_);
    }

    Finite<P> pow(int n) const;

    Finite<P> rev();

};

template<unsigned int P>
Finite<P> &Finite<P>::operator+=(const Finite<P> &a) {
    value_ += a.value_;
    value_ = normalize(value_);
    return *this;
}

template<unsigned int P>
Finite<P> &Finite<P>::operator-=(const Finite &a) {
    value_ -= a.value_;
    value_ = normalize(value_);
    return *this;
}

template<unsigned int P>
Finite<P> &Finite<P>::operator*=(const Finite &a) {
    long long res = static_cast<long long>(value_) *  static_cast<long long>(a.value_);
    value_ = normalize(res);
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

template<int X>
struct checking_ {
};

template<>
struct checking_<0> {
    void fake_function() {
    }
};

template<int N, int M>
struct get_equal_ {
    get_equal_() {
        checking_<N - M> x;
        x.fake_function();
    }


};

template<size_t N, size_t M, typename Field = Rational>
class Matrix {
private:
    template<size_t SZ, typename Field_>
    class Matrix_row {

        //Matrix_row<SZ, Field_> &operator=(Matrix_row<SZ, Field_>) {}
    private:
        std::vector<Field_> v_;

    public:
        Matrix_row() {
            v_.resize(SZ);
        }

        explicit Matrix_row(std::vector<Field_> &g) {
            v_ = g;
        }

        Field_ &operator[](size_t x) {
            return v_[x];
        }

        const Field_ &operator[](size_t x) const {
            return v_[x];
        }
    };


    void universal_gauss_();

    std::vector<Matrix_row<M, Field>> v_;
    Field det_;
    size_t rank_ = -1;
    typedef std::vector<std::vector<Field>> subm_;
    typedef std::vector<std::vector<subm_>> parts_;

    static subm_ copy_part_(const subm_ &from, size_t xs, size_t xe, size_t ys, size_t ye) {
        size_t n = from.size();
        subm_ to;
        to.resize(n / 2);
        for (size_t i = 0; xs + i < xe; ++i) {
            to[i].resize(n / 2);
            for (size_t j = 0; ys + j < ye; ++j) {
                to[i][j] = from[xs + i][ys + j];
            }
        }
        return to;
    }

    static subm_ strassen_substract_(const subm_ &a, const subm_ &b) {
        size_t n = a.size();
        subm_ ans;
        ans.resize(n);
        for (size_t i = 0; i < n; ++i) {
            ans[i].resize(n);
            for (size_t j = 0; j < n; ++j) {
                ans[i][j] = a[i][j] - b[i][j];
            }
        }
        return ans;
    }

    static subm_ strassen_add_(const subm_ &a, const subm_ &b) {
        size_t n = a.size();
        subm_ ans;
        ans.resize(n);
        for (size_t i = 0; i < n; ++i) {
            ans[i].resize(n);
            for (size_t j = 0; j < n; ++j) {
                ans[i][j] = a[i][j] + b[i][j];
            }
        }
        return ans;
    }

    static subm_ strassen_multiply_(const subm_ &a, const subm_ &b) {
        size_t n = a.size();
        subm_ c;
        if (n <= 32) {
            c.resize(n);
            for (size_t i = 0; i < n; i++) {
                c[i].resize(n);
                for (size_t j = 0; j < n; j++) {
                    for (size_t k = 0; k < n; k++) {
                        c[i][j] += (a[i][k] * b[k][j]);
                    }
                }
            }
            return c;
        } else {
            parts_ m_a;
            parts_ m_b;
            m_a.resize(2);
            m_b.resize(2);
            for (size_t i = 0; i < 2; ++i) {
                m_a[i].resize(2);
                m_b[i].resize(2);
                for (size_t j = 0; j < 2; ++j) {
                    m_a[i][j].resize(n / 2);
                    m_b[i][j].resize(n / 2);
                    for (size_t k = 0; k < n / 2; ++k) {
                        m_a[i][j][k].resize(n / 2);
                        m_b[i][j][k].resize(n / 2);
                        for (size_t l = 0; l < n / 2; ++l) {
                            m_a[i][j][k][l] = a[i * (n / 2) + k][j * (n / 2) + l];
                            m_b[i][j][k][l] = b[i * (n / 2) + k][j * (n / 2) + l];
                        }
                    }
                }
            }
            subm_ m1 = strassen_multiply_(strassen_add_(m_a[0][0], m_a[0][1]), strassen_add_(m_b[0][0], m_b[0][1]));
            subm_ m2 = strassen_multiply_(strassen_add_(m_a[1][0], m_a[1][1]), m_b[0][0]);
            subm_ m3 = strassen_multiply_(m_a[0][0], strassen_substract_(m_b[0][1], m_b[1][1]));
            subm_ m4 = strassen_multiply_(m_a[1][1], strassen_substract_(m_b[1][0], m_b[0][0]));
            subm_ m5 = strassen_multiply_(strassen_add_(m_a[0][0], m_a[0][1]), m_b[1][1]);
            subm_ m6 = strassen_multiply_(strassen_substract_(m_a[1][0], m_a[0][0]),
                                          strassen_add_(m_b[0][0], m_b[0][1]));
            subm_ m7 = strassen_multiply_(strassen_substract_(m_a[0][1], m_a[1][1]),
                                          strassen_add_(m_b[1][0], m_b[1][1]));
            parts_ m_c;
            m_c.resize(2, std::vector<subm_>(2));
            m_c[0][0] = strassen_substract_(strassen_add_(m1, m4), strassen_add_(m5, m7));
            m_c[0][1] = strassen_add_(m3, m5);
            m_c[1][0] = strassen_add_(m2, m4);
            m_c[1][1] = strassen_add_(strassen_substract_(m1, m2), strassen_add_(m3, m6));
            subm_ ans;
            ans.resize(n);
            for (size_t i = 0; i < n; ++i) {
                ans[i].resize(n);
                for (size_t j = 0; j < n; ++j) {
                    ans[i][j] = m_c[i / (n / 2)][j / (n / 2)][i % (n / 2)][j % (n / 2)];
                }
            }
            return ans;
        }
    }

public:
    Matrix<N, M, Field>(const std::initializer_list<std::initializer_list<int>> &arr) {
        v_.resize(N);
        size_t i = 0;
        for (auto it : arr) {
            size_t j = 0;

            for (auto item : it) {
                (*this)[i][j] = Field(item);
                ++j;
            }
            ++i;
        }
    }

    explicit Matrix(const std::vector<std::vector<Field>> &x);

    Matrix() {
        v_.resize(N);
        for (size_t i = 0; i < std::min(N, M); i++) {
            (*this)[i][i] = Field(1);
        }
    }

    Matrix<N, M, Field> &operator+=(const Matrix<N, M, Field> &a);

    Matrix<N, M, Field> operator+(const Matrix<N, M, Field> &a) const;

    Matrix<N, M, Field> &operator-=(const Matrix<N, M, Field> &a);

    Matrix<N, M, Field> operator-(const Matrix<N, M, Field> &a) const;

    bool operator!=(const Matrix<N, M, Field> &a) const;

    /*template<size_t K>
    friend Matrix<N, K, Field> operator*=(Matrix<N, M, Field> &a, const Matrix<M, K, Field> &b){
        get_equal_<N,M> x;
        get_equal_<M,K> y;
        *this = a * b;
        return *this;
    }*/

    template<size_t K>
    friend Matrix<N, K, Field> operator*(const Matrix<N, M, Field> &a, const Matrix<M, K, Field> &b) {
        size_t max_zs = std::max(N, std::max(M, K));
        size_t pow = 2;
        while (pow < max_zs) {
            pow *= 2;
        }
        subm_ aa;
        subm_ bb;
        subm_ cc;
        aa.resize(pow, std::vector<Field>(pow));
        bb.resize(pow, std::vector<Field>(pow));
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                aa[i][j] = a[i][j];
            }
        }
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < K; ++j) {
                bb[i][j] = b[i][j];
            }
        }
        cc = strassen_multiply_(aa, bb);
        Matrix<N, K, Field> ans;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < K; ++j) {
                ans[i][j] = cc[i][j];
            }
        }

        return ans;
    }

    Field det() const;

    Matrix<M, N, Field> transposed() const;

    size_t rank() const;

    Matrix<N, M, Field> &invert();

    Matrix<N, M, Field> inverted() const;

    Field trace() const;

    std::vector<Field> getRow(size_t x) const;

    std::vector<Field> getColumn(size_t y) const;

    Matrix_row<M, Field> &operator[](size_t x) {
        return v_[x];
    };

    const Matrix_row<M, Field> &operator[](size_t x) const {
        return v_[x];
    };

    Matrix<N, M, Field> &operator*=(const Field &x);

    Matrix<N, M, Field> &operator/=(const Field &x);

    Matrix<N, M, Field> operator*(const Field &x) const;

    Matrix<N, M, Field> operator/(const Field &x) const;
};

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field>::Matrix(const std::vector<std::vector<Field>> &x) {
    v_ = x;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> &Matrix<N, M, Field>::operator+=(const Matrix<N, M, Field> &a) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            v_[i][j] += a[i][j];
        }
    }
    return *this;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::operator+(const Matrix<N, M, Field> &a) const {
    Matrix<N, M, Field> res = *this;
    res += a;
    return res;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> &Matrix<N, M, Field>::operator-=(const Matrix<N, M, Field> &a) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            v_[i][j] -= a[i][j];
        }
    }
    return *this;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::operator-(const Matrix<N, M, Field> &a) const {
    Matrix<N, M, Field> res = *this;
    res -= a;
    return res;
}

template<size_t N, size_t M, typename Field>
Matrix<M, N, Field> Matrix<N, M, Field>::transposed() const {
    Matrix<M, N, Field> res;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            res[j][i] = (*this)[i][j];
        }
    }
    return res;
}

template<size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::universal_gauss_() {
    size_t rank = 0;
    Field det = 1;
    for (size_t i = 0; i < M; ++i) {
        size_t k = rank;
        for (size_t j = rank + 1; j < N; ++j)
            if ((((*this))[j][i]) > (((*this))[k][i]))
                k = j;
        if ((((*this))[k][i]) == Field(0)) {
            det = Field(0);
            continue;
            //здесь тут
        } else {
            //rank++;
        }
        Matrix_row<M, Field> x = (*this)[rank];
        (*this)[rank] = (*this)[k];
        (*this)[k] = x;
        //std::swap((*this)[i], (*this)[k]);
        if (i != k) {
            det = -det;
        }
        det *= (*this)[rank][i];
        for (size_t j = i + 1; j < M; ++j) {
            (*this)[rank][j] /= (*this)[rank][i];
        }
        (*this)[rank][i] /= (*this)[rank][i];
        for (size_t j = 0; j < N; ++j) {
            if (j != rank && ((*this)[j][i]) != Field(0)) {
                for (size_t k = i + 1; k < M; ++k) {
                    (*this)[j][k] -= (*this)[rank][k] * (*this)[j][i];
                }
                //(*this)[j][i] = Field(0);
            }
        }
        rank++;

    }

    det_ = det;
    rank_ = rank;

}

template<size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::det() const {
    Matrix<N, M, Field> res = *this;
    //res = res.transposed();
    Field det = Field(1);
    for (size_t i = 0; i < N; ++i) {
        size_t k = i;
        for (size_t j = i + 1; j < N; ++j)
            if ((res[j][i]) > (res[k][i]))
                k = j;
        if ((res[k][i]) == Field(0)) {
            //det = Field(0);
            //break;
        }
        std::swap(res[i], res[k]);
        if (i != k)
            det = -det;
        det *= res[i][i];
        for (size_t j = i + 1; j < N; ++j)
            res[i][j] /= res[i][i];
        for (size_t j = 0; j < N; ++j)
            if (j != i /*&& (res[j][i]) != Field(0)*/)
                for (size_t k = i + 1; k < N; ++k)
                    res[j][k] -= res[i][k] * res[j][i];
    }
    return det;
}

template<size_t N, size_t M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getRow(size_t x) const {
    std::vector<Field> ans;
    for (size_t i = 0; i < M; i++) {
        ans.push_back((*this)[x][i]);
    }
    return ans;
}

template<size_t N, size_t M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getColumn(size_t y) const {
    std::vector<Field> ans;
    for (size_t i = 0; i < N; i++) {
        ans.push_back((*this)[i][y]);
    }
    return ans;
}

template<size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::trace() const {
    get_equal_<N, M> x;
    Field tr;
    for (size_t i = 0; i < N; i++) {
        tr += (*this)[i][i];
    }
    return tr;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> &Matrix<N, M, Field>::invert() {
    get_equal_<N, M> test;
    Matrix<N, M, Field> one;
    for (size_t i = 0; i < N; i++) {
        one[i][i] = Field(1);
    }

    for (size_t i = 0; i < M; ++i) {
        size_t k = i;
        for (size_t j = i + 1; j < N; ++j) {
            if ((((*this))[j][i]) > (((*this))[k][i])) {
                k = j;
            }
        }
        if (((*this)[k][i]) == Field(0)) {
            break;
            //throw exception
        }
        //Matrix_row<M, Field> x = (*this)[i];
        //(*this)[i] = (*this)[k];
        // (*this)[k] = x;
        std::swap((*this)[i], (*this)[k]);
        std::swap(one[i], one[k]);
        for (size_t j = 0; j < N; ++j) {
            if (j > i) {
                (*this)[i][j] /= (*this)[i][i];
            }
            one[i][j] /= (*this)[i][i];
        }
        (*this)[i][i] /= (*this)[i][i];
        for (size_t j = 0; j < N; ++j) {
            if (j != i && ((*this)[j][i]) != Field(0)) {
                for (size_t k = 0; k < M; ++k) {
                    if (i != k) {
                        (*this)[j][k] -= ((*this)[i][k] * (*this)[j][i]);
                    }
                    //std::cout << i << " " << j << " " << k << std::endl;

                    one[j][k] -= (one[i][k] * (*this)[j][i]);
                }
                //(*this)[j][i] = Field(0);

            }
        }
        /*for (int i = 0; i < N; i++) {
             for (int j = 0; j < M; j++) {
                 std::cout << ((*this)[i][j]).value_ << " ";
             }
             std::cout << std::endl;
         }
         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;*/
    }

    *this = one;
    return *this;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::inverted() const {
    Matrix<N, M, Field> tmp = (*this);
    return tmp.invert();
}

template<size_t N, size_t M, typename Field>
size_t Matrix<N, M, Field>::rank() const {
    Matrix<N, M, Field> x = (*this);
    x.universal_gauss_();
    return x.rank_;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Field &x, const Matrix<N, M, Field> &matrix) {
    return matrix * x;
}

template<size_t N, typename Field>
Matrix<N, N, Field> &operator*=(Matrix<N, N, Field> &first, const Matrix<N, N, Field> &second) {
    return first = (first * second);
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> &Matrix<N, M, Field>::operator*=(const Field &x) {
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            (*this)[i][j] *= x;
        }
    }
    return *this;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::operator*(const Field &x) const {
    Matrix<N, M, Field> ans = *this;
    return ans *= x;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> &Matrix<N, M, Field>::operator/=(const Field &x) {
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            *(this)[i][j] /= x;
        }
    }
    return *this;
}

template<size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::operator/(const Field &x) const {
    Matrix<N, M, Field> ans = *this;
    return ans /= x;

}

template<size_t N, size_t M, typename Field>
bool Matrix<N, M, Field>::operator!=(const Matrix<N, M, Field> &a) const {
    bool fl = false;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            if ((*this)[i][j] != a[i][j]) {
                fl = true;
            }
        }
    }
    return fl;
}

template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;
#endif //MATRIX_ZP_MATRIX_H
