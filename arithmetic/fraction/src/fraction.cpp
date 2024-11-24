#include "../include/fraction.h"

fraction::fraction(
    big_integer &&numerator,
    big_integer &&denominator):
        _numerator(std::forward<big_integer>(numerator)),
        _denominator(std::forward<big_integer>(denominator))
{
    if(_denominator.is_equal_to_zero())
        throw std::logic_error("Division by zero");

    big_integer gcd = big_integer::euclidean_algorithm(this->_numerator, this->_denominator);
    this->_numerator /= gcd;
    this->_denominator /= gcd;

}

fraction::fraction(const big_integer &numerator, const big_integer &denominator): _numerator(numerator), _denominator(denominator)
{

}

fraction::~fraction()
{

}

fraction::fraction(
    fraction const &other):
        _numerator(other._numerator),
        _denominator(other._denominator)
{

}

fraction &fraction::operator=(
    fraction const &other)
{
    if(*this != other)
    {
        this->_numerator = other._numerator;
        this->_denominator = other._denominator;
    }
    return *this;
}

fraction::fraction(
    fraction &&other) noexcept:
        _numerator(std::move(other._numerator)),
        _denominator(std::move(other._denominator))
{

}

fraction &fraction::operator=(
    fraction &&other) noexcept
{
    if(*this != other)
    {
        this->_numerator = std::move(other._numerator);
        this->_denominator = std::move(other._denominator);
    }
    return *this;
}

fraction &fraction::operator+=(
    fraction const &other)
{
    this->_numerator = this->_numerator * other._denominator + other._numerator * this->_denominator;
    this->_denominator *= other._denominator;

    if(this->_numerator.sign() == -1 && this->_denominator.sign() == -1
       || this->_numerator.sign() == 1 && this->_denominator.sign() == -1)
    {
        this->_numerator = std::move(-this->_numerator);
        this->_denominator = std::move(-this->_denominator);
    }


    big_integer gcd = big_integer::euclidean_algorithm(this->_numerator, this->_denominator);
    this->_numerator /= gcd;
    this->_denominator /= gcd;

    return *this;
}

fraction fraction::operator+(
    fraction const &other) const
{
    return fraction(*this) += other;
}

fraction &fraction::operator-=(
    fraction const &other)
{
    this->_numerator = this->_numerator * other._denominator - other._numerator * this->_denominator;
    this->_denominator *= other._denominator;

    if(this->_numerator.sign() == -1 && this->_denominator.sign() == -1
       || this->_numerator.sign() == 1 && this->_denominator.sign() == -1)
    {
        this->_numerator = std::move(-this->_numerator);
        this->_denominator = std::move(-this->_denominator);
    }

    big_integer gcd = big_integer::euclidean_algorithm(this->_numerator, this->_denominator);
    this->_numerator /= gcd;
    this->_denominator /= gcd;

    return *this;
}

fraction fraction::operator-(
    fraction const &other) const
{
    return fraction(*this) -= other;
}

fraction &fraction::operator*=(
    fraction const &other)
{
//    std::cout<<this->_numerator * other._numerator<<std::endl;
//    std::cout<<this->_denominator * other._denominator<<std::endl;
    this->_numerator *= other._numerator;
    this->_denominator *= other._denominator;

    if(this->_numerator.sign() == -1 && this->_denominator.sign() == -1
       || this->_numerator.sign() == 1 && this->_denominator.sign() == -1)
    {
        this->_numerator = std::move(-this->_numerator);
        this->_denominator = std::move(-this->_denominator);
    }

    big_integer gcd = big_integer::euclidean_algorithm(this->_numerator, this->_denominator);
    this->_numerator /= gcd;
    this->_denominator /= gcd;

    return *this;
}

fraction fraction::operator*(
    fraction const &other) const
{
    return fraction(*this) *= other;
}

fraction &fraction::operator/=(
    fraction const &other)
{
    this->_numerator *= other._denominator;
    this->_denominator *= other._numerator;

    if(this->_numerator.sign() == -1 && this->_denominator.sign() == -1
    || this->_numerator.sign() == 1 && this->_denominator.sign() == -1)
    {
        this->_numerator = std::move(-this->_numerator);
        this->_denominator = std::move(-this->_denominator);
    }

    big_integer gcd = big_integer::euclidean_algorithm(this->_numerator, this->_denominator);
    this->_numerator /= gcd;
    this->_denominator /= gcd;

    return *this;
}

fraction fraction::operator/(
    fraction const &other) const
{
    return fraction(*this) /= other;
}

bool fraction::operator==(
    fraction const &other) const
{
    if(this->_numerator == other._numerator && this->_denominator == other._denominator)
        return true;

    return false;
}

bool fraction::operator!=(
    fraction const &other) const
{
    return !(*this == other);
}

bool fraction::operator>=(
    fraction const &other) const
{
    fraction copied_this(*this);
    fraction copied_other(other);

    copied_this._denominator *= other._denominator;
    copied_this._numerator *= other._denominator;

    copied_other._denominator *= this->_denominator;
    copied_other._numerator *= this->_denominator;

    return copied_this._numerator >= copied_other._numerator;
}

bool fraction::operator>(
    fraction const &other) const
{
    fraction copied_this(*this);
    fraction copied_other(other);

    copied_this._denominator *= other._denominator;
    copied_this._numerator *= other._denominator;

    copied_other._denominator *= this->_denominator;
    copied_other._numerator *= this->_denominator;

    return copied_this._numerator > copied_other._numerator;
}

bool fraction::operator<=(
    fraction const &other) const
{
    fraction copied_this(*this);
    fraction copied_other(other);

    copied_this._denominator *= other._denominator;
    copied_this._numerator *= other._denominator;

    copied_other._denominator *= this->_denominator;
    copied_other._numerator *= this->_denominator;

    return copied_this._numerator <= copied_other._numerator;
}

bool fraction::operator<(
    fraction const &other) const
{
    fraction copied_this(*this);
    fraction copied_other(other);

    copied_this._denominator *= other._denominator;
    copied_this._numerator *= other._denominator;

    copied_other._denominator *= this->_denominator;
    copied_other._numerator *= this->_denominator;

    return copied_this._numerator < copied_other._numerator;
}

std::ostream &operator<<(
    std::ostream &stream,
    fraction const &obj)
{
    if(obj._numerator.sign() == 1)
        stream<<'+';
    stream<<obj._numerator<<"/"<<obj._denominator;
    return stream;
}

std::istream &operator>>(
    std::istream &stream,
    fraction &obj)
{
    char sign;
    char division_symbol;
    std::string numerator;
    std::string denominator;

    stream>>sign>>numerator>>division_symbol>>denominator;

    if(sign == '-')
        numerator = sign + numerator;

    obj = fraction(big_integer(numerator), big_integer(denominator));

    return stream;
}

fraction fraction::factorial(
        size_t number)
{
    fraction result(big_integer("1"), big_integer("1"));
    for(size_t iter = 2; iter <= number; ++iter)
    {
        big_integer factorial_element(std::to_string(iter));
        result._numerator *= factorial_element;
    }
    return result;
}

fraction fraction::abs(
        const fraction number)
{
    big_integer new_numerator = number._numerator;
    if(new_numerator.sign() == -1)
        new_numerator = -new_numerator;

    return fraction(new_numerator, number._denominator);
}

fraction fraction::sin(
    fraction const &epsilon) const
{

    /*fraction copied_this(*this);
    fraction pi_2(big_integer("44"), big_integer("7"));
    size_t iter = 0;

    while(copied_this > pi_2)
    {
        copied_this /= pi_2;
        iter++;
    }
    if(iter!= 0)
        copied_this = *this / pi_2 / fraction(big_integer(std::to_string(iter)), big_integer("1"));*/

    auto sin_raw = [](size_t n, const fraction &x) -> fraction {
        fraction term = x.pow(2 * n + 1) / factorial(2 * n + 1);
        if (n & 1) {
            term._numerator = std::move(-term._numerator);
        }
        return term;
    };

    fraction sum(big_integer("0"), big_integer("1"));
    fraction prev_sum = sum;
    fraction current(big_integer("0"), big_integer("1"));
    size_t n = 0;
    do {
        std::cout<<n<<std::endl;
        prev_sum = sum;
        current = sin_raw(n, *this);
        sum += current;
        ++n;
    } while (abs(current) >= epsilon && n <= 200);


    return sum;
}
fraction fraction::cos(
    fraction const &epsilon) const
{
    auto cos_raw = [](size_t n, fraction current) -> fraction
    {
        fraction counted = current.pow(2 * n) / factorial(2 * n);

        if(n & 1)
               counted._numerator = std::move(-counted._numerator);

        return counted;
    };

    fraction sum(big_integer("0"), big_integer("1"));
    fraction prev_sum(sum);
    fraction current(big_integer("0"), big_integer("1"));
    size_t n = 0;
    do
    {

        prev_sum = sum;
        current = cos_raw(n, *this);
        sum += current;
        ++n;
    }while(abs(current) >= epsilon && n <= 200);

    return sum;
}

fraction fraction::tg(
    fraction const &epsilon) const
{
    return this->sin(epsilon) / this->cos(epsilon);
}

fraction fraction::ctg(
    fraction const &epsilon) const
{
    return this->cos(epsilon) / this->sin(epsilon);
}

fraction fraction::sec(
    fraction const &epsilon) const
{
    fraction result = this->cos(epsilon);

    big_integer tmp = std::move(result._denominator);
    result._denominator = std::move(result._numerator);
    result._numerator = std::move(tmp);

    if(result._numerator.sign() == -1 && result._denominator.sign() == -1
       || result._numerator.sign() == 1 && result._denominator.sign() == -1)
    {
        result._numerator = std::move(-result._numerator);
        result._denominator = std::move(-result._denominator);
    }

    return result;
}

fraction fraction::cosec(
    fraction const &epsilon) const
{
    fraction result = this->sin(epsilon);

    big_integer tmp = std::move(result._denominator);
    result._denominator = std::move(result._numerator);
    result._numerator = std::move(tmp);


    if(result._numerator.sign() == -1 && result._denominator.sign() == -1
       || result._numerator.sign() == 1 && result._denominator.sign() == -1)
    {
        result._numerator = std::move(-result._numerator);
        result._denominator = std::move(-result._denominator);
    }

    return result;
}

fraction fraction::arcsin(
    fraction const &epsilon) const
{
    if(abs(*this) > fraction(big_integer("1"), big_integer("1")))
        throw std::logic_error("The raw converges at -1 <= x <= 1");

    auto arcsin_raw = [](size_t n, fraction current) -> fraction
    {
        big_integer numerator1("1");
        big_integer denominator("1");

        for(size_t i = 2; i <= 2 * n + 1; ++i)
        {
            big_integer multiplier(std::to_string(i));

            if(i & 1)
                numerator1 *= multiplier;
            else
                denominator *= multiplier;
        }

        fraction multiplier1(std::move(numerator1), std::move(denominator));

        fraction multiplier2 = current.pow(2* n + 1) / fraction(big_integer(std::to_string(2 * n +1)), big_integer("1")).pow(2);

        return multiplier1 * multiplier2;
    };

    fraction prev_summ = arcsin_raw(0,*this);

    fraction current = arcsin_raw(1, *this);
    fraction result = prev_summ + current;

    size_t iter = 2;
    while(abs(current) > epsilon && iter <= 200)
    {
        prev_summ = result;
        current = arcsin_raw(iter, *this);
        result += current;
        iter++;
    }
    return result;
}

fraction fraction::arccos(
    fraction const &epsilon) const
{
    if(abs(*this) > fraction(big_integer("1"), big_integer("1")))
        throw std::logic_error("The raw converges at -1 <= x <= 1");

    return fraction(big_integer("1"), big_integer("1")).arcsin(epsilon) - this->arcsin(epsilon);
}

fraction fraction::arctg(
    fraction const &epsilon) const
{

    auto arctg_raw = [](size_t n, fraction current) -> fraction
    {
        fraction result = current.pow(2 * n + 1) / fraction(big_integer(std::to_string(2 * n + 1)), big_integer("1"));
        if( n & 1)
            result._numerator = std::move(-result._numerator);

        return result;
    };

    fraction prev_summ = arctg_raw(0, *this);

    fraction current = arctg_raw(1, *this);
    fraction result = prev_summ + current;
    size_t iter = 2;

    while(abs(current) >= epsilon && iter <= 30)
    {
        prev_summ = result;
        current = arctg_raw(iter, *this);
        result += current;
        iter++;
    }
    return result;
}

fraction fraction::arcctg(
    fraction const &epsilon) const
{

    return this->arctg(epsilon)
    + fraction(big_integer("2"), big_integer("1"))
    * fraction(big_integer("1"), big_integer("1")).arctg(epsilon);
}

fraction fraction::arcsec(
    fraction const &epsilon) const
{
    big_integer new_numerator("0");
    big_integer new_denominator("0");

    if(this->_numerator.sign() == -1)
    {
        new_numerator = -this->_denominator;
        new_denominator = -this->_numerator;
    }
    else
    {
        new_numerator = this->_denominator;
        new_denominator = this->_numerator;
    }

    return fraction(new_numerator,new_denominator).arccos(epsilon);
}

fraction fraction::arccosec(
    fraction const &epsilon) const
{
    big_integer new_numerator("0");
    big_integer new_denominator("0");

    if(this->_numerator.sign() == -1)
    {
        new_numerator = -this->_denominator;
        new_denominator = -this->_numerator;
    }
    else
    {
        new_numerator = this->_denominator;
        new_denominator = this->_numerator;
    }

    return fraction(new_numerator,new_denominator).arcsin(epsilon);
}

fraction fraction::pow(
    size_t degree) const
{

    auto binary_pow = [](fraction base, size_t degree)->fraction
    {
        fraction result(big_integer("1"), big_integer("1"));

        while(degree)
        {
            if(degree & 1)
                result *= base;
            base *= base;
            degree >>= 1;
        }
        return result;
    };

    return binary_pow(*this, degree);
}

fraction fraction::root(
    size_t degree,
    fraction const &epsilon) const
{
    if (this->_numerator.sign() == -1 && (degree & 1) == 0) {
        throw std::domain_error("Cannot compute the even root of a negative fraction.");
    }


    fraction copy = *this;
    bool swapped = false;
    if (copy._numerator > copy._denominator)
    {
        std::swap(copy._numerator, copy._denominator);
        swapped = true;
    }

    fraction alpha = fraction(big_integer("1"), big_integer(std::to_string(degree)));
    copy -= fraction(big_integer("1"), big_integer("1"));

    fraction sum(big_integer("1"), big_integer("1"));
    fraction prev = fraction(big_integer("2"), big_integer("1")) * epsilon;
    int iteration = 1;

    fraction precompute = alpha;
    while (abs(prev) > epsilon && iteration <= 100) {
        std::cout<<iteration<<std::endl;
        prev = precompute;
        prev *= copy.pow(iteration);
        prev /= factorial(iteration);
        sum += prev;

        precompute *= (alpha - fraction(big_integer(std::to_string(iteration)), big_integer("1")));
        ++iteration;
    }

    if (swapped) {
        std::swap(sum._denominator, sum._numerator);
    }

    return sum;
}

fraction fraction::log2(
    fraction const &epsilon) const
{
    if(this->_numerator.sign() == -1)
        throw std::logic_error("You cannot calculate logarithm of negative argument");
    if(this->_numerator == this->_denominator)
        return fraction(big_integer("0"), big_integer("1"));

    return ln(epsilon)/fraction(big_integer("2"), big_integer("1")).ln(epsilon);
}

fraction fraction::ln(
    fraction const &epsilon) const
{
    if(this->_numerator.sign() == -1)
        throw std::logic_error("You cannot calculate logarithm of negative argument");

    if(this->_numerator == this->_denominator)
        return fraction(big_integer("0"), big_integer("1"));


    fraction x = *this;
    if (x < fraction(big_integer("1"), big_integer("1")))
    {
        big_integer tmp = std::move(x._numerator);
        x._numerator = std::move(x._denominator);
        x._denominator = std::move(tmp);
    }
    fraction z = (fraction(big_integer("1"), big_integer("1")) - x) / x;

    auto ln_raw = [](size_t n, fraction z) -> fraction
    {
        fraction counted = z.pow(n)
                / fraction(big_integer(std::to_string(n)), big_integer("1"));

        if((n - 1) & 1)
            counted._numerator = std::move(-counted._numerator);

        return counted;
    };

    fraction sum(big_integer("0"), big_integer("1"));
    fraction current = fraction(big_integer("1"), big_integer("1"));
    size_t n = 1;
    do {
//        std::cout<<sum<<" "<<current<<std::endl;
        current = ln_raw(n, z);
        sum += current;
        n++;
    } while (abs(current) >= epsilon && n <= 200);

    if(*this > fraction(big_integer("1"), big_integer("1")))
        sum._numerator = std::move(-sum._numerator);
    return sum;

}

fraction fraction::lg(
    fraction const &epsilon) const
{
    if(this->_numerator.sign() == -1)
        throw std::logic_error("You cannot calculate logarithm of negative argument");
    if(this->_numerator == this->_denominator)
        return fraction(big_integer("0"), big_integer("1"));

    return ln(epsilon) / fraction(big_integer("10"), big_integer("1")).ln(epsilon);
}