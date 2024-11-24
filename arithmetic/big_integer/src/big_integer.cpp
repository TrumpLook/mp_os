#include <cstring>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstdint>

#include "big_integer.h"

void big_integer::clear()
{
    _oldest_digit = 0;
//    delete[] _other_digits;
    deallocate_with_guard(_other_digits);
    _other_digits = nullptr;
}

void big_integer::copy_from(
        big_integer const &other)
{
    _oldest_digit = other._oldest_digit;
    _other_digits = nullptr;
    if (other._other_digits == nullptr)
    {
        return;
    }

    _other_digits = reinterpret_cast<unsigned int*>(allocate_with_guard(sizeof(unsigned int), *other._other_digits));

//    _other_digits = new unsigned int[*other._other_digits];
    std::memcpy(_other_digits, other._other_digits, sizeof(unsigned int) * (*other._other_digits));
}

void big_integer::initialize_from(
        int const *digits,
        size_t digits_count)
{
    if (digits == nullptr)
    {
        throw std::logic_error("pointer to digits array must not be nullptr");
    }

    if (digits_count == 0)
    {
        throw std::logic_error("digits array length must  be GT 0");
    }

    _oldest_digit = digits[digits_count - 1];
//    _other_digits = (digits_count == 1
//                     ? nullptr
//                     : new unsigned int[digits_count]);

    _other_digits = (digits_count == 1
            ? nullptr
            : reinterpret_cast<unsigned int*>(allocate_with_guard(sizeof(unsigned int), *_other_digits)));

    if (_other_digits == nullptr)
    {
        return;
    }

    *_other_digits = (unsigned int)digits_count;

    std::memcpy(_other_digits + 1, digits, sizeof(unsigned int) * (digits_count - 1));
}

void big_integer::initialize_from(
        std::vector<int> const &digits,
        size_t digits_count)
{
    _other_digits = nullptr;

    if (digits.empty() || digits_count == 0)
    {
        throw std::logic_error("std::vector<int> of digits should not be empty");
    }

    _oldest_digit = digits[digits_count - 1];

    if (digits_count == 1)
    {
        return;
    }

//    _other_digits = new unsigned int[digits_count];

    _other_digits = reinterpret_cast<unsigned int*>(allocate_with_guard(sizeof(unsigned int), digits_count));
    *_other_digits = digits_count;

    for (auto i = 0; i < digits_count - 1; ++i)
    {
        _other_digits[1 + i] = *reinterpret_cast<unsigned int const *>(&digits[i]);
    }
}

void big_integer::initialize_from(
        const std::vector<unsigned int> &digits,
        size_t digits_count)
{
    _other_digits = nullptr;

    if(digits.empty() || digits_count == 0)
    {
        throw std::logic_error("std::vector<unsigned int> of digits shouldn't be empty");
    }

    _oldest_digit = digits[digits_count - 1];

    if(digits_count == 1)
        return;

//    _other_digits = new unsigned int[digits_count];

    _other_digits = reinterpret_cast<unsigned int*>(allocate_with_guard(sizeof(unsigned int), digits_count));
    *_other_digits = digits_count;

    for(int i = 0; i < digits_count - 1; ++i)
        _other_digits[1 + i] = digits[i];

}

void big_integer::initialize_from(
        std::string const &value,
        size_t base)
{
    std::vector<unsigned int> converted = convert_to_base(value, base);
    initialize_from(converted, converted.size());
}

void big_integer::print_byte(
        std::ostream &stream,
        unsigned char byte_value)
{
    for (int i = 0; i < 8; i++)
    {
        stream << ((byte_value >> (7 - i)) & 1);
    }
}

void big_integer::dump_int_value(
        std::ostream &stream,
        int value)
{
    auto *p = (unsigned char *)&value;
    for (int i = 0; i < sizeof(int); i++)
    {
        print_byte(stream, *p++);
        stream << ' ';
    }
}

big_integer &big_integer::change_sign()
{
    _oldest_digit ^= (1 << ((sizeof(int) << 3) - 1));

    return *this;
}

std::vector<unsigned int> big_integer::convert_to_base(std::string const &biiiiiiiiiiig_number, size_t base)
{
    bool is_negative = false;

    int pos = 0;

    if(biiiiiiiiiiig_number[0] == '-')
    {
        is_negative = true;
        pos = 1;
    }
    auto divide_str_on_int = [](std::string str, int position = 0) -> std::vector<unsigned int>
    {
        std::vector<unsigned int> result;
        size_t converted = 0;

        unsigned int max_int = -1;
        size_t base = static_cast<size_t>(max_int) + 1;

        while(position != str.length())
        {
            std::string next_number_to_divide("");

            while(converted < base)
            {
                if(position == str.length()) break;

                converted = converted * 10 + (str[position] - '0');
                position ++;
            }
            if(position == str.length())
            {
                if(converted >= base)
                {
                    result.push_back(converted % base);
                    converted /= base;
                }
                result.push_back(converted);
                return  result;
            }

            while(position != str.length())
            {

                if(converted >= base)
                {
                    next_number_to_divide.push_back(converted / base + '0');
                    converted %= base;
                }
                else next_number_to_divide.push_back('0');

                if(position != str.length()) converted = converted * 10 + (str[position] - '0');

                position++;
            }
            if(converted >= base)
            {
                next_number_to_divide.push_back(converted / base + '0');
                converted %= base;
            }

            else next_number_to_divide.push_back('0');

            result.push_back(converted);
            str = std::move(next_number_to_divide);
            converted = 0;
            position = 0;
        }

        return result;
    };

    std::vector<unsigned int> result = divide_str_on_int(biiiiiiiiiiig_number, pos);

    if((result[result.size() - 1] & (1 << ((sizeof(unsigned int) << 3) - 1))) != 0)
    {
        result.push_back(0);
    }


    if(result.size() == 1 && result.back() == 0)
        return result;

    if(is_negative)
    {
        result[result.size() - 1] |= 1<<(sizeof(int) * 8 - 1);
    }
    return result;

}

inline int big_integer::get_digits_count() const noexcept
{
    return static_cast<int>(_other_digits == nullptr
                            ? 1
                            : *_other_digits);
}

inline int big_integer::sign() const noexcept
{
    if (is_equal_to_zero())
    {
        return 0;
    }

    return 1 - (static_cast<int>((*reinterpret_cast<unsigned int const *>(&_oldest_digit) >> ((sizeof(int) << 3) - 1))) << 1);
}

inline bool big_integer::is_equal_to_zero() const noexcept
{
    return _oldest_digit == 0 && _other_digits == nullptr;
}

inline unsigned int big_integer::get_digit(
        int position) const noexcept
{
    if (_other_digits == nullptr)
    {
        return position == 0
               ? _oldest_digit
               : 0;
    }

    int const digits_count = get_digits_count();

    if (position < digits_count - 1)
    {
        return _other_digits[position + 1];
    }

    if (position == digits_count - 1)
    {
        return _oldest_digit;
    }

    return 0;
}

allocator *big_integer::get_allocator() const noexcept
{
    return this->_allocator;
}
big_integer::big_integer(std::string const& digits, allocator* allocator) : _allocator(allocator)
{
    initialize_from(digits, digits.size());
}

big_integer::big_integer(
        int const *digits,
        size_t digits_count,
        allocator* allocator) : _allocator(allocator)
{
    initialize_from(digits, digits_count);
}

big_integer::big_integer(
        std::vector<int> const &digits,
        allocator* allocator) : _allocator(allocator)
{
    initialize_from(digits, digits.size());
}

big_integer::big_integer(
        std::string const &value,
        size_t base,
        allocator* allocator) : _allocator(allocator)
{
    initialize_from(value, base);
}

big_integer::big_integer(
        big_integer const &other)
{
    this->_allocator = other._allocator;
    copy_from(other);
}

big_integer &big_integer::operator=(
        big_integer const &other)
{
    this->_allocator = other._allocator;
    if (this != &other)
    {
        clear();
        copy_from(other);
    }

    return *this;
}

big_integer::~big_integer()
{
    clear();
}

void big_integer::remove_leading_zeros()
{
    if(_other_digits == nullptr)
        return;
    if(_oldest_digit == 0)
    {
        std::vector<unsigned int> new_digits;
        for(int i = 0; i <= get_digits_count() - 1; ++i)
            new_digits.push_back(get_digit(i));

        while(new_digits.back() == 0)
            new_digits.pop_back();

        new_digits.push_back(0);
/*        for(int i = 0; i < new_digits.size(); ++i)
            std::cout<<new_digits[i]<<" ";
        std::cout<<std::endl;*/

        clear();
        initialize_from(new_digits, new_digits.size());
    }
}

big_integer big_integer::euclidean_algorithm(big_integer first_number, big_integer second_number)
{
    big_integer tmp("0");

    while (second_number != big_integer("0"))
    {
        tmp = first_number % second_number;
        first_number = second_number;
        second_number = tmp;
        second_number.remove_leading_zeros();
/*        for(int i = 0; i < second_number.get_digits_count() - 1;  ++i)
            std::cout<<second_number.get_digit(i)<<" ";
        std::cout<<std::endl;
        std::cout<<first_number<<" "<<second_number<<std::endl<<std::endl;*/

    }
    return first_number;
}

big_integer &big_integer::operator+=(
        big_integer const &other)
{
    if (other.is_equal_to_zero())
    {
        return *this;
    }

    if (is_equal_to_zero())
    {
        return *this = other;
    }

    if (sign() == -1)
    {
        change_sign();
        *this += -other;
        return change_sign();
    }

    if (other.sign() == -1)
    {
        return *this -= -other;
    }

    auto const first_value_digits_count = get_digits_count();
    auto const second_value_digits_count = other.get_digits_count();
    auto const digits_count = std::max(first_value_digits_count, second_value_digits_count);

    unsigned int operation_result = 0;

    constexpr int shift = sizeof(unsigned int) << 2;
    constexpr int mask = (1 << shift) - 1;

    std::vector<unsigned int> result_digits(digits_count + 1);

    for (int i = 0; i < digits_count; ++i)
    {
        unsigned int first_value_digit = get_digit(i);
        unsigned int second_value_digit = other.get_digit(i);
        result_digits[i] = 0;

        for (int j = 0; j < 2; ++j)
        {
            operation_result += (first_value_digit & mask) + (second_value_digit & mask);
            first_value_digit >>= shift;
            second_value_digit >>= shift;
            result_digits[i] |= ((operation_result & mask) << shift * j);
            operation_result >>= shift;
        }
    }

    while(result_digits.back() == 0)
        result_digits.pop_back();

   if(operation_result)
        result_digits.back() += operation_result;

    if(result_digits[result_digits.size() - 1] & (1 << (sizeof(unsigned int) * 8 - 1)))
        result_digits.push_back(0);


    auto result_digits_count = result_digits.size();

    clear();
    initialize_from(result_digits, result_digits_count);

    return *this;
}

big_integer big_integer::operator+(
        big_integer const &other) const
{
    return big_integer(*this) += other;
}

big_integer &big_integer::operator-=(
        big_integer const &other)
{
    if(this->is_equal_to_zero() && !other.is_equal_to_zero())
    {
        big_integer tmp(other);

        tmp.change_sign();
        *this = std::move(tmp);
         return *this;
    }

    else if(!this->is_equal_to_zero() && other.is_equal_to_zero())
        return *this;

    else if(this->sign() == -1 && other.sign() == 1)
    {
        this->change_sign();
        *this += other;
        this->change_sign();
        return *this;
    }
    else if(this->sign() == 1 && other.sign() == -1)
    {
        big_integer tmp(other);
        tmp.change_sign();

        *this += tmp;
        return *this;
    }

    else if(this->sign() == -1 && other.sign() == -1)
    {
        big_integer tmp(other);

        tmp.change_sign();

        tmp += *this;

        *this = std::move(tmp);
        return *this;
    }

    else if(this->sign() == 1 && other.sign() == 1)
    {
        if(*this < other)
        {
           return *this = -(other - *this);
        }
        else
        {
             bool need_to_borrow = false;

            std::vector<unsigned int> result;
            for(int i = 0; i < this->get_digits_count(); ++i)
            {
                auto number_one = this->get_digit(i);
                auto number_two = i < other.get_digits_count() ? other.get_digit(i) : 0;

                unsigned int operation_result = number_one - number_two - need_to_borrow;

                need_to_borrow = number_one < number_two;

                result.push_back(operation_result);
            }
            big_integer tmp("0");
            while(result[result.size() - 1] == 0)
                result.pop_back();

            if(result.empty())
                result.push_back(0);

            if((result[result.size() - 1] & (1 << ((sizeof(unsigned int) << 3) - 1))) != 0)
                result.push_back(0);

            tmp.initialize_from(result, result.size());
            *this = std::move(tmp);

            return *this;
        }
    }
    else
    {
        *this = std::move(big_integer("0"));
        return *this;
    }

}

big_integer big_integer::operator-(
        big_integer const &other) const
{
    return big_integer(*this) -= other;
}

big_integer big_integer::operator-() const
{
    return big_integer(*this).change_sign();
}

big_integer &big_integer::operator*=(
        big_integer const &other)
{
    big_integer copied_other(other);

    bool _is_result_negative;

    if(this->sign() == -1 && copied_other.sign() == 1 || this->sign() == 1 && copied_other.sign() == -1)
    {
        _is_result_negative = true;
        this->sign() == -1 ? this->change_sign() : copied_other.change_sign();
    }
    else if(this->sign() == -1 && copied_other.sign() == -1)
    {
        _is_result_negative = false;
        this->change_sign();
        copied_other.change_sign();
    }
    else _is_result_negative = false;



    int const size_of_this = this->get_digits_count();
    int const size_of_other = copied_other.get_digits_count();

    constexpr int shift = sizeof(unsigned int) << 2;
    constexpr unsigned int mask = (1U << shift) - 1;

    big_integer result("0");

    for(int i = 0; i < 2 * size_of_this; ++i)
    {

        unsigned int remainder = 0;



        unsigned int first_number_half;

        if(i % 2 == 0)
        {
            auto number = this->get_digit(i / 2);
            first_number_half = number & mask;
        }
        else
        {
            auto number = this->get_digit(i / 2);
            first_number_half = (number >> shift) & mask;
        }

        for(int j = 0; j < 2 * size_of_other; ++j)
        {
            std::vector<unsigned int> digits_array;
            unsigned int second_number_half;
            if(j % 2 == 0)
            {
                auto number = copied_other.get_digit(j / 2);
                second_number_half = number & mask;
            }
            else
            {
                auto number = copied_other.get_digit(j / 2);
                second_number_half = (number >> shift) & mask;
            }

            unsigned int operation_result = (first_number_half * second_number_half + remainder) & mask;
            remainder = (first_number_half * second_number_half + remainder) >> shift;


            digits_array.push_back(operation_result);

            big_integer multiply_result("0");
            multiply_result.initialize_from(digits_array, digits_array.size());

            multiply_result <<= (shift * (i + j));
            result += multiply_result;


        }
        if(remainder)
        {
            std::vector<unsigned int> remainder_vector(1);
            remainder_vector[0] = remainder;
            big_integer add_remainder("0");
            add_remainder.initialize_from(remainder_vector, remainder_vector.size());



            add_remainder <<= (shift * (2 * size_of_other + i));
            result += add_remainder;


        }
    }

    *this = std::move(result);

    if(_is_result_negative) this->change_sign();
    return *this;
}

big_integer big_integer::operator*(
        big_integer const &other) const
{
    return big_integer(*this) *= other;
}

big_integer &big_integer::operator/=(
        big_integer const &other)
{
    if(other.is_equal_to_zero())
        throw std::logic_error("NOPE, DUDE, YOU CANNOT DIVIDE ON ZERO");

    else if(this->is_equal_to_zero())
        return *this = std::move(big_integer("0"));

    std::stringstream this_str;
    this_str<<*this;
    std::stringstream other_str;
    other_str<<other;


    /**this = std::move(big_integer(this_str.str()));*/

    bool is_result_negative = false;
    big_integer copied_other(other);

    if(this->sign() ^ other.sign())
    {
        if(this->sign() == 1)
            copied_other.change_sign();
        else
            this->change_sign();
        is_result_negative = true;
    }

    else
    {
        if(this->sign() == -1)
        {
            copied_other.change_sign();
            this->change_sign();
        }
    }

    if(this->_oldest_digit == 1 && this->_other_digits == nullptr || copied_other._oldest_digit == 1 && copied_other._other_digits ==
                                                                                                        nullptr)
    {
        if(copied_other._oldest_digit == 1 && copied_other._other_digits == nullptr)
        {
            if(is_result_negative)
                this->change_sign();

            return *this;
        }
    }

    {
        big_integer tmp(other);
        big_integer one("1");
        int shift = 0;

        while((tmp & one) != one)
        {
            tmp>>=1;
            shift++;
        }

        if(tmp == one)
        {
            *this >>= shift;
            if(is_result_negative)
                this->change_sign();
            return *this;
        }
    }

    auto compare_vectors = [](std::vector<unsigned int>& vector_1, std::vector<unsigned int>& vector_2) -> int
    {
        if(!vector_1.empty())
        {
            if(vector_1.front() == 0)
            {
                std::vector<unsigned int> new_vec;
                int last_zero = 0;
                for(int i = 1; i < vector_1.size(); ++i)
                {
                    if(vector_1[i] == 0)
                        last_zero++;
                    else
                        break;
                }

                for(int i = last_zero + 1; i < vector_1.size(); ++i)
                    new_vec.push_back(vector_1[i]);

                vector_1 = std::move(new_vec);
            }
        }

        if(!vector_2.empty())
        {
            if(vector_2.front() == 0)
            {
                std::vector<unsigned int> new_vec;
                int last_zero = 0;
                for(int i = 1; i < vector_2.size(); ++i)
                {
                    if(vector_2[i] == 0)
                        last_zero++;
                    else
                        break;
                }

                for(int i = last_zero + 1; i < vector_2.size(); ++i)
                    new_vec.push_back(vector_2[i]);

                vector_2 = std::move(new_vec);
            }
        }


        if(vector_1.size() > vector_2.size())
            return 1;

        else if(vector_1.size() < vector_2.size())
            return -1;

        else
        {
            for(size_t i = 0; i < vector_1.size(); ++i)
            {
                if(vector_1[i] == vector_2[i])
                    continue;

                else if(vector_1[i] > vector_2[i])
                    return 1;

                else
                    return -1;
            }
        }
        return 0;
    };

    auto multiply_vector_on_int = [](unsigned int number, std::vector<unsigned int> const & big_number) -> std::vector<unsigned int>
    {
        std::vector<unsigned int> copied_number(big_number);
        std::vector<unsigned int> result;

        unsigned int const max_int = -1;
        size_t base = static_cast<size_t>(max_int) + 1;

        std::reverse(copied_number.begin(), copied_number.end());

        unsigned int remainder = 0;

        for(int i = 0; i < copied_number.size(); ++i)
        {
            size_t multiplication_result = static_cast<size_t>(copied_number[i]) * static_cast<size_t>(number) + remainder;

            result.push_back(multiplication_result % base);

            remainder = multiplication_result / base;
        }

        if(remainder)
            result.push_back(remainder);

        std::reverse(result.begin(), result.end());

        return result;
    };

    auto swap_or_add_bits = [] (unsigned int& number, bool swap, unsigned int& left_number, unsigned int& right_number)
    {
        if(left_number == right_number)
            number = left_number;
        else if(swap)
        {
            size_t mid = static_cast<size_t>(left_number) + static_cast<size_t>(right_number);
            mid /= 2;
            right_number = number;
            number = mid;
        }
        else
        {
            size_t mid = static_cast<size_t>(left_number) + static_cast<size_t>(right_number);
            mid /= 2;
            left_number = number;
            number = mid + 1;
        }

    };

    std::vector<unsigned int> result;
    std::vector<unsigned int> dividend;


    std::vector<unsigned int> divider;
    for(int i = copied_other.get_digits_count() - 1; i >= 0; --i)
        if(!(divider.empty() && copied_other.get_digit(i) == 0))
            divider.push_back(copied_other.get_digit(i));


    for(int i = this->get_digits_count() - 1; i >=0; --i)
    {
        if(compare_vectors(dividend, divider) == -1)
        {

            if(!(dividend.empty() && this->get_digit(i) == 0))
                dividend.push_back(this->get_digit(i));

            result.push_back(0);
            if(compare_vectors(dividend, divider) == -1) continue;
        }

        unsigned int divide_result = 1 << ((sizeof(unsigned int) << 3) - 1);
        unsigned int left = 1;
        unsigned int right = -1;


        while(true)
        {
            std::vector<unsigned int> multiplication_result;

            multiplication_result = multiply_vector_on_int(divide_result, divider);

            auto comparasion_result = compare_vectors(multiplication_result, dividend);

            if(comparasion_result == 1)
            {
                swap_or_add_bits(divide_result, true, left, right);
            }

            else if(comparasion_result == -1)
            {
                big_integer big_dividend("0");
                std::vector<unsigned int> copied_dividend_vector(dividend);
                std::reverse(copied_dividend_vector.begin(), copied_dividend_vector.end());

                if((copied_dividend_vector[copied_dividend_vector.size() - 1] & (1 << ((sizeof(unsigned int) << 3) - 1))) != 0)
                    copied_dividend_vector.push_back(0);

                big_dividend.initialize_from(copied_dividend_vector, copied_dividend_vector.size());

                std::reverse(multiplication_result.begin(), multiplication_result.end());
                big_integer big_mult_result("0");


                if((multiplication_result[multiplication_result.size() - 1] & (1 << ((sizeof(unsigned int) << 3) - 1))) != 0)
                    multiplication_result.push_back(0);


                big_mult_result.initialize_from(multiplication_result, multiplication_result.size());

                if(multiplication_result.back() == 0)
                    multiplication_result.pop_back();


                big_integer remainder = big_dividend - big_mult_result;
//                std::cout<<remainder<<" "<<big_dividend<<" "<<copied_other<<std::endl;


                if(remainder >= copied_other)
                {
                    swap_or_add_bits(divide_result, false, left, right);
                }
                else
                {
                    result.push_back(static_cast<int>(divide_result));
                    dividend.clear();
                    if(!remainder.is_equal_to_zero())
                    {
                        for(int j = remainder.get_digits_count() - 1; j >= 0; --j)
                            dividend.push_back(remainder.get_digit(j));
                    }

                    break;
                }
            }

            else
            {
                result.push_back(static_cast<int>(divide_result));
                std::reverse(dividend.begin(), dividend.end());
                std::reverse(multiplication_result.begin() ,multiplication_result.end());

                big_integer number_1("0");
                number_1.initialize_from(dividend, dividend.size());

                big_integer number_2("0");
                number_2.initialize_from(multiplication_result, multiplication_result.size());

                number_1 -= number_2;

                dividend.clear();

                if(!number_1.is_equal_to_zero())
                {
                    for(int j = number_1.get_digits_count() - 1; j >= 0; --j)
                        dividend.push_back(number_1.get_digit(j));
                }
                break;
            }
        }

        dividend.push_back(this->get_digit(i - 1));
    }

    std::reverse(result.begin(), result.end());

    big_integer result_big("0");
    result_big.initialize_from(result, result.size());

    if(is_result_negative)
        result_big.change_sign();

    *this = std::move(result_big);
    return *this;
}

big_integer big_integer::operator/(
        big_integer const &other) const
{
    return big_integer(*this) /= other;
}

big_integer &big_integer::operator%=(
        big_integer const &other)
{
    if(other.is_equal_to_zero())
        throw std::logic_error("NOPE, DUDE, YOU CANNOT DIVIDE ON ZERO");

    else if(this->is_equal_to_zero())
        return *this = std::move(big_integer("0"));

    bool is_result_negative = false;
    big_integer copied_other(other);

    if(this->sign() ^ other.sign())
    {
        if(this->sign() == 1)
            copied_other.change_sign();
        else
            this->change_sign();
        is_result_negative = true;
    }

    else
    {
        if(this->sign() == -1)
        {
            copied_other.change_sign();
            this->change_sign();
        }
    }

    if(*this < copied_other)
        return *this;

    if(copied_other._oldest_digit == 1 && copied_other._other_digits == nullptr)
    {
        if(copied_other._oldest_digit == 1 && copied_other._other_digits == nullptr)
        {
            if(is_result_negative)
                this->change_sign();

            *this = std::move(big_integer("0"));

            return *this;
        }
    }

    auto compare_vectors = [](std::vector<unsigned int> & vector_1, std::vector<unsigned int> & vector_2) -> int
    {
        if(!vector_1.empty())
        {
            if(vector_1.front() == 0)
            {
                std::vector<unsigned int> new_vec;
                int last_zero = 0;
                for(int i = 1; i < vector_1.size(); ++i)
                {
                    if(vector_1[i] == 0)
                        last_zero++;
                    else
                        break;
                }

                for(int i = last_zero + 1; i < vector_1.size(); ++i)
                    new_vec.push_back(vector_1[i]);

                vector_1 = std::move(new_vec);
            }
        }

        if(!vector_2.empty())
        {
            if(vector_2.front() == 0)
            {
                std::vector<unsigned int> new_vec;
                int last_zero = 0;
                for(int i = 1; i < vector_2.size(); ++i)
                {
                    if(vector_2[i] == 0)
                        last_zero++;
                    else
                        break;
                }

                for(int i = last_zero + 1; i < vector_2.size(); ++i)
                    new_vec.push_back(vector_2[i]);

                vector_2 = std::move(new_vec);
            }
        }

        if(vector_1.size() > vector_2.size())
            return 1;

        else if(vector_1.size() < vector_2.size())
            return -1;

        else
        {
            for(size_t i = 0; i < vector_1.size(); ++i)
            {
                if(vector_1[i] == vector_2[i])
                    continue;

                else if(vector_1[i] > vector_2[i])
                    return 1;

                else
                    return -1;
            }
        }
        return 0;
    };

    auto multiply_vector_on_int = [](unsigned int number, std::vector<unsigned int> const & big_number) -> std::vector<unsigned int>
    {
        std::vector<unsigned int> copied_number(big_number);
        std::vector<unsigned int> result;

        unsigned int const max_int = -1;
        size_t base = static_cast<size_t>(max_int) + 1;

        std::reverse(copied_number.begin(), copied_number.end());

        unsigned int remainder = 0;

        for(int i = 0; i < copied_number.size(); ++i)
        {
            size_t multiplication_result = static_cast<size_t>(copied_number[i]) * static_cast<size_t>(number) + remainder;

            result.push_back(multiplication_result % base);

            remainder = multiplication_result / base;
        }

        if(remainder)
            result.push_back(remainder);

        std::reverse(result.begin(), result.end());

        return result;
    };

    auto swap_or_add_bits = [] (unsigned int& number, bool swap, unsigned int& left_number, unsigned int& right_number)
    {
        if(left_number == right_number)
            number = left_number;
        else if(swap)
        {
            size_t mid = static_cast<size_t>(left_number) + static_cast<size_t>(right_number);
            mid /= 2;
            right_number = number;
            number = mid;
        }
        else
        {
            size_t mid = static_cast<size_t>(left_number) + static_cast<size_t>(right_number);
            mid /= 2;
            left_number = number;
            number = mid + 1;
        }

    };

    std::vector<unsigned int> result;
    std::vector<unsigned int> dividend;

    big_integer remainder(*this);

    std::vector<unsigned int> divider;
    for(int i = copied_other.get_digits_count() - 1; i >= 0; --i)
        if(!(divider.empty() && copied_other.get_digit(i) == 0))
            divider.push_back(copied_other.get_digit(i));


    for(int i = this->get_digits_count() - 1; i >=0; --i)
    {
        if(compare_vectors(dividend, divider) == -1)
        {
            if(!(dividend.empty() && this->get_digit(i) == 0))
                dividend.push_back(this->get_digit(i));

            result.push_back(0);
            if(compare_vectors(dividend, divider) == -1) continue;
        }

        unsigned int divide_result = 1 << ((sizeof(unsigned int) << 3) - 1);
        unsigned int left = 1;
        unsigned int right = -1;


        while(true)
        {
            std::vector<unsigned int> multiplication_result;

            multiplication_result = multiply_vector_on_int(divide_result, divider);

            auto comparasion_result = compare_vectors(multiplication_result, dividend);

            if(comparasion_result == 1)
            {
                swap_or_add_bits(divide_result, true, left, right);
            }

            else if(comparasion_result == -1)
            {
                big_integer big_dividend("0");
                std::vector<unsigned int> copied_dividend_vector(dividend);
                std::reverse(copied_dividend_vector.begin(), copied_dividend_vector.end());

                if((copied_dividend_vector[copied_dividend_vector.size() - 1] & (1 << ((sizeof(unsigned int) << 3) - 1))) != 0)
                    copied_dividend_vector.push_back(0);

                big_dividend.initialize_from(copied_dividend_vector, copied_dividend_vector.size());

                std::reverse(multiplication_result.begin(), multiplication_result.end());
                big_integer big_mult_result("0");


                if((multiplication_result[multiplication_result.size() - 1] & (1 << ((sizeof(unsigned int) << 3) - 1))) != 0)
                    multiplication_result.push_back(0);


                big_mult_result.initialize_from(multiplication_result, multiplication_result.size());

                if(multiplication_result.back() == 0)
                    multiplication_result.pop_back();


                remainder = big_dividend - big_mult_result;

                if(remainder >= copied_other)
                {
                    swap_or_add_bits(divide_result, false, left, right);
                }
                else
                {
                    result.push_back(static_cast<int>(divide_result));
                    dividend.clear();
                    if(!remainder.is_equal_to_zero())
                    {
                        for(int j = remainder.get_digits_count() - 1; j >= 0; --j)
                            dividend.push_back(remainder.get_digit(j));
                    }

                    break;
                }
            }

            else
            {
                remainder = big_integer("0");

                result.push_back(static_cast<int>(divide_result));
                std::reverse(dividend.begin(), dividend.end());
                std::reverse(multiplication_result.begin() ,multiplication_result.end());

                big_integer number_1("0");
                number_1.initialize_from(dividend, dividend.size());

                big_integer number_2("0");
                number_2.initialize_from(multiplication_result, multiplication_result.size());

                number_1 -= number_2;

                dividend.clear();

                if(!number_1.is_equal_to_zero())
                {
                    for(int j = number_1.get_digits_count() - 1; j >= 0; --j)
                        dividend.push_back(number_1.get_digit(j));
                }
                break;
            }
        }

        dividend.push_back(this->get_digit(i - 1));
    }

    *this = std::move(remainder);
    return *this;
}

big_integer big_integer::operator%(
        big_integer const &other) const
{
    return big_integer(*this) %= other;
}

bool big_integer::operator==(
        big_integer const &other) const
{

    if(this->sign() != other.sign()) return false;

    auto this_size = this->get_digits_count();
    auto other_size = other.get_digits_count();

    if(this_size != other_size) return false;

    for(int i = this_size -1; i >= 0; --i)
        if(this->get_digit(i) != other.get_digit(i)) return false;

    return true;
}

bool big_integer::operator!=(
        big_integer const &other) const
{
    return !(*this == other);
}

bool big_integer::operator<(
        big_integer const &other) const
{
    if(this->sign() != other.sign())
        return this->sign() == -1;

    else
    {
        big_integer copied_other(other);
        big_integer copied_this(*this);

        if(other.sign() == -1)
        {
            copied_this.change_sign();
            copied_other.change_sign();
        }

        auto this_digits_count = this->get_digits_count();
        auto other_digits_count = other.get_digits_count();

        if(this_digits_count > other_digits_count)
            return this->sign() == -1;

        else if(this_digits_count < other_digits_count)
            return this->sign() == 1;


        else
        {


            for(int i = this_digits_count - 1; i >= 0; --i)
            {
                if(copied_this.get_digit(i) == copied_other.get_digit(i))
                    continue;

                if(copied_this.get_digit(i) > copied_other.get_digit(i))
                    return this->sign() == -1;

                else if(copied_this.get_digit(i) < copied_other.get_digit(i))
                    return this->sign() == 1;
            }
            return false;
        }

    }
}

bool big_integer::operator<=(
        big_integer const &other) const
{
    if(this->sign() != other.sign())
        return this->sign() == -1;

    else
    {
        big_integer copied_other(other);
        big_integer copied_this(*this);

        if(other.sign() == -1)
        {
            copied_this.change_sign();
            copied_other.change_sign();
        }

        auto this_digits_count = this->get_digits_count();
        auto other_digits_count = other.get_digits_count();

        if(this_digits_count > other_digits_count)
            return this->sign() == -1;

        else if(this_digits_count < other_digits_count)
            return this->sign() == 1;


        else
        {


            for(int i = this_digits_count - 1; i >= 0; --i)
            {
                if(copied_this.get_digit(i) == copied_other.get_digit(i))
                    continue;

                if(copied_this.get_digit(i) > copied_other.get_digit(i))
                    return this->sign() == -1;

                else if(copied_this.get_digit(i) < copied_other.get_digit(i))
                    return this->sign() == 1;
            }
            return true;
        }

    }
}

bool big_integer::operator>(
        big_integer const &other) const
{
    if(this->sign() != other.sign())
        return this->sign() == 1;

    else
    {
        big_integer copied_other(other);
        big_integer copied_this(*this);

        if(other.sign() == -1)
        {
            copied_this.change_sign();
            copied_other.change_sign();
        }

        auto this_digits_count = this->get_digits_count();
        auto other_digits_count = other.get_digits_count();

        if(this_digits_count > other_digits_count)
            return this->sign() == 1;

        else if(this_digits_count < other_digits_count)
            return this->sign() == -1;


        else
        {


            for(int i = this_digits_count - 1; i >= 0; --i)
            {
                if(copied_this.get_digit(i) == copied_other.get_digit(i))
                    continue;

                if(copied_this.get_digit(i) > copied_other.get_digit(i))
                    return this->sign() == 1;

                else if(copied_this.get_digit(i) < copied_other.get_digit(i))
                    return this->sign() == -1;
            }
            return false;
        }

    }
}

bool big_integer::operator>=(
        big_integer const &other) const
{
    if(this->sign() != other.sign())
        return this->sign() == 1;

    else
    {
        big_integer copied_other(other);
        big_integer copied_this(*this);

        if(other.sign() == -1)
        {
            copied_this.change_sign();
            copied_other.change_sign();
        }

        auto this_digits_count = this->get_digits_count();
        auto other_digits_count = other.get_digits_count();

        if(this_digits_count > other_digits_count)
            return this->sign() == 1;

        else if(this_digits_count < other_digits_count)
            return this->sign() == -1;


        else
        {
            for(int i = this_digits_count - 1; i >= 0; --i)
            {
                if(copied_this.get_digit(i) == copied_other.get_digit(i))
                    continue;

                if(copied_this.get_digit(i) > copied_other.get_digit(i))
                    return this->sign() == 1;

                else if(copied_this.get_digit(i) < copied_other.get_digit(i))
                    return this->sign() == -1;
            }
            return true;
        }

    }
}

big_integer big_integer::operator~() const
{
    auto digits = get_digits_count();


//    int* new_digits = new int[digits + 1];

    int* new_digits = reinterpret_cast<int*>(allocate_with_guard(sizeof(int), digits + 1));
    new_digits[0] = digits;

    for(int i = 1; i < digits; ++i)
    {
        auto current = get_digit(i);
        current = ~current;
        new_digits[i] = *reinterpret_cast<int*>(&current);
    }

    big_integer result(new_digits, *new_digits);
    return result;

}

big_integer &big_integer::operator&=(
        big_integer const &other)
{
    int const size_1 = this->get_digits_count();
    int const size_2 = other.get_digits_count();

    int const new_size = std::min(size_1, size_2);

    std::vector<unsigned int> new_digits;


    for(int i = 0; i < new_size; ++i)
    {
        new_digits.push_back(this->get_digit(i) & other.get_digit(i));
    }

    this->clear();
    this->initialize_from(new_digits, new_size);

    return *this;
}

big_integer big_integer::operator&(
        big_integer const &other) const
{
    return big_integer(*this) &= other;
}

big_integer &big_integer::operator|=(
        big_integer const &other)
{
    int const size_1 = this->get_digits_count();
    int const size_2 = other.get_digits_count();

    int const new_size =  std::max(size_1, size_2);
//    int* new_digits = new int[new_size + 1];

    int* new_digits = reinterpret_cast<int*>(allocate_with_guard(sizeof(int), new_size + 1));

    *new_digits = new_size;

    for(int i = 0; i < new_size; ++i)
    {
        if(i != size_1 && i != size_2) new_digits[i] = this->get_digit(i) | other.get_digit(i);
        else
        {
             new_digits[i] = i == size_1 ? other.get_digit(i) : this->get_digit(i);
        }
    }
    this->clear();
    this->initialize_from(new_digits, new_size);

    return *this;
}

big_integer big_integer::operator|(
        big_integer const &other) const
{
    return big_integer(*this) |= other;
}

big_integer &big_integer::operator^=(
        big_integer const &other)
{
    int const size_1 = this->get_digits_count();
    int const size_2 = other.get_digits_count();

    int const new_size =  std::max(size_1, size_2);
//    int* new_digits = new int[new_size + 1];

    int* new_digits = reinterpret_cast<int*>(allocate_with_guard(sizeof(int), new_size + 1));

    *new_digits = new_size;

    for(int i = 0; i < new_size; ++i)
    {
        if(i != size_1 && i != size_2) new_digits[i] = this->get_digit(i) ^ other.get_digit(i);
        else
        {
            new_digits[i] = i == size_1 ? other.get_digit(i) : this->get_digit(i);
        }
    }
    this->clear();
    this->initialize_from(new_digits, new_size);

    return *this;
}

big_integer big_integer::operator^(
        big_integer const &other) const
{
    return big_integer(*this) ^= other;
}

big_integer &big_integer::operator<<=(
        size_t shift_value)
{
    if (is_equal_to_zero() || shift_value == 0)
    {
        return *this;
    }

    auto value_sign = sign();
    if (value_sign == -1)
    {
        change_sign();
    }

    auto const added_by_shift_at_other_digits_digits_count = shift_value / (sizeof(unsigned int) << 3);
    shift_value %= (sizeof(unsigned int) << 3);

    auto added_by_shift_at_oldest_digit_digits_count = 0;
    if (_oldest_digit != 0)
    {
        unsigned int oldest_digit = *reinterpret_cast<unsigned int *>(&_oldest_digit);
        int oldest_value_bit_index = 0;
        while (oldest_digit != 1)
        {
            oldest_digit >>= 1;
            ++oldest_value_bit_index;
        }

        if (oldest_value_bit_index + shift_value > (sizeof(int) << 3) - 2)
        {
            ++added_by_shift_at_oldest_digit_digits_count;
        }
    }

    if (added_by_shift_at_oldest_digit_digits_count != 0 || added_by_shift_at_other_digits_digits_count != 0)
    {
        auto const added_digits_count = added_by_shift_at_oldest_digit_digits_count + added_by_shift_at_other_digits_digits_count;

        if (_other_digits == nullptr)
        {
//            _other_digits = new unsigned int[added_digits_count + 1];
            _other_digits = reinterpret_cast<unsigned int*>(allocate_with_guard(sizeof(unsigned int), added_digits_count  + 1));
            *_other_digits = added_digits_count + 1;
            std::memset(_other_digits + 1, 0, sizeof(unsigned int) * (added_digits_count - 1));
            if (added_by_shift_at_oldest_digit_digits_count != 0)
            {
                _other_digits[*_other_digits - 1] = _oldest_digit;
                _oldest_digit = 0;
            }
            else
            {
                _other_digits[*_other_digits - 1] = 0;
            }
        }
        else
        {
//            auto *new_digits = new unsigned int[added_digits_count + *_other_digits];
            unsigned int* new_digits = reinterpret_cast<unsigned int*>(allocate_with_guard(sizeof(unsigned int), added_digits_count + *_other_digits));
            std::memset(new_digits + 1, 0, sizeof(unsigned int) * added_digits_count);
            if (added_by_shift_at_oldest_digit_digits_count != 0)
            {
                new_digits[added_digits_count + *_other_digits - 1] = _oldest_digit;
                _oldest_digit = 0;
            }
            std::memcpy(new_digits + 1 + added_by_shift_at_other_digits_digits_count, _other_digits + 1, sizeof(unsigned int) * (*_other_digits - 1));
            *new_digits = *_other_digits + added_digits_count;

//            delete[] _other_digits;
            deallocate_with_guard(_other_digits);
            _other_digits = new_digits;
        }
    }

    if (shift_value != 0)
    {
        auto const digits_count = get_digits_count();
        unsigned int part_to_move_to_next_digit = 0;
        for (auto i = 0; i < digits_count; ++i)
        {
            auto digit_value = get_digit(i);
            auto *digit_address = i == digits_count - 1
                                  ? reinterpret_cast<unsigned int *>(&_oldest_digit)
                                  : _other_digits + 1 + i;
            *digit_address <<= shift_value;
            *digit_address |= part_to_move_to_next_digit;
            part_to_move_to_next_digit = digit_value >> ((sizeof(unsigned int) << 3) - shift_value);
        }
    }

    if (value_sign == -1)
    {
        change_sign();
    }

    return *this;
}

big_integer big_integer::operator<<(
        size_t shift_value) const
{
    return big_integer(*this) <<= shift_value;
}

big_integer &big_integer::operator>>=(
        size_t shift_value)
{
    if (is_equal_to_zero() || shift_value == 0) {
        return *this;
    }

    auto value_sign = sign();
    if (value_sign == -1) {
        change_sign();
    }

    auto const remove_digits_count = shift_value / (sizeof(unsigned int) << 3);
    shift_value %= (sizeof(unsigned int) << 3);

    if (remove_digits_count > 0) {
        if (remove_digits_count >= *_other_digits) {
//            delete[] _other_digits;
            deallocate_with_guard(_other_digits);
            _other_digits = nullptr;
            _oldest_digit = 0;
            return *this;
        }

        auto new_size = *_other_digits - remove_digits_count;
//        unsigned int *new_digits = new unsigned int[new_size + 1];
        unsigned int *new_digits = reinterpret_cast<unsigned int*>(allocate_with_guard(sizeof(unsigned int), new_size + 1));
        std::memcpy(new_digits + 1, _other_digits + 1 + remove_digits_count, sizeof(unsigned int) * (new_size - 1));
        *new_digits = new_size;
//        delete[] _other_digits;
        deallocate_with_guard(_other_digits);
        _other_digits = new_digits;
    }

    if (shift_value != 0) {
        unsigned int part_to_move_to_previous_digit = 0;
        auto const digits_count = get_digits_count();
        for (auto i = digits_count - 1; i >= 0; --i) {
            auto *digit_address = i == digits_count - 1
                                  ? reinterpret_cast<unsigned int *>(&_oldest_digit)
                                  : _other_digits + 1 + i;
            unsigned int current_digit = get_digit(i);
            *(digit_address) >>= shift_value;
            *(digit_address) |= part_to_move_to_previous_digit;
            part_to_move_to_previous_digit = (current_digit << ((sizeof(unsigned int) << 3) - shift_value));
        }
    }

    if (value_sign == -1) {
        change_sign();
    }

    return *this;
}

big_integer big_integer::operator>>(
        size_t shift_value) const
{
    return big_integer(*this) >>= shift_value;
}

std::ostream &operator<<(
        std::ostream &stream,
        big_integer const &value)
{

    unsigned int max_int = -1;
    size_t base = static_cast<size_t>(max_int) + 1;
    std::string base_in_str = std::to_string(base);

    big_integer integer(value);

    if(value.sign() == -1)
        integer.change_sign();



    auto multiply_strings = [](const std::string &num1, const std::string &num2) {
        int n1 = num1.size();
        int n2 = num2.size();
        std::vector<int> result(n1 + n2, 0);


        std::string num1Rev(num1.rbegin(), num1.rend());
        std::string num2Rev(num2.rbegin(), num2.rend());


        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                result[i + j] += (num1Rev[i] - '0') * (num2Rev[j] - '0');
                result[i + j + 1] += result[i + j] / 10;
                result[i + j] %= 10;
            }
        }


        int i = result.size() - 1;
        while (i > 0 && result[i] == 0)
            --i;


        std::string res;
        while (i >= 0)
            res.push_back(result[i--] + '0');

        return res;
    };

    auto sum_two_numbers_in_string = [](std::string& result, std::string& str_to_add)->std::string
    {
        std::reverse(result.begin(), result.end());
        std::reverse(str_to_add.begin(), str_to_add.end());

        std::string tmp("");

        const int size = std::max(result.length(), str_to_add.length());

        const int size_diff = result.length() - str_to_add.length();

        if(size_diff != 0)
            size_diff < 0 ? result.resize(result.length()-size_diff, '0') : str_to_add.resize(str_to_add.length() + size_diff, '0');

        int next_degree = 0;

        for(int i = 0; i < size; ++i)
        {
            int first_digit = str_to_add[i] - '0';
            int second_digit = result[i] - '0';

            int sum_result = first_digit + second_digit + next_degree;

            tmp.push_back((sum_result % 10) + '0');
            next_degree = sum_result / 10;
        }
        if(next_degree) tmp.push_back('1');

        std::reverse(tmp.begin(), tmp.end());
        return tmp;
    };

    std::string converted("");

    int const digits_count = integer.get_digits_count();

    for(int i = 0; i < digits_count - 1; ++i)
    {
        auto number =  integer.get_digit(i);

        std::string current_digit_in_str("1");

        for(int j = 0; j < i; ++j)
            current_digit_in_str = multiply_strings(current_digit_in_str, base_in_str);


        current_digit_in_str = multiply_strings(current_digit_in_str, std::to_string(number));

        std::string tmp = sum_two_numbers_in_string(converted, current_digit_in_str);
        converted = std::move(tmp);
    }

    {
        auto last_digit = integer._oldest_digit;

        std::string current_digit_in_str("1");

        for(int i = 0; i < digits_count - 1; ++i)
            current_digit_in_str = multiply_strings(current_digit_in_str, base_in_str);

        current_digit_in_str = multiply_strings(current_digit_in_str, std::to_string(last_digit));

        std::string tmp = sum_two_numbers_in_string(converted, current_digit_in_str);
        converted = std::move(tmp);
    }


    if(value.sign() == -1) converted = "-" + converted;
    stream<<(converted);
    return stream;
}

std::istream &operator>>(
        std::istream &stream,
        big_integer &value)
{
    std::string source;
    stream >> source;
    big_integer tmp(source);

    value = std::move(tmp);

    return stream;

}

big_integer &big_integer::trivial_multiplication::multiply(big_integer &first_multiplier,
                                                           const big_integer &second_multiplier) const
                                                           {
                                                                return first_multiplier*=second_multiplier;
                                                           }

big_integer &big_integer::trivial_division::divide(big_integer &dividend, const big_integer &divisor,
                                                   big_integer::multiplication_rule multiplication_rule) const
                                                   {
                                                        return dividend /= divisor;
                                                   }

big_integer &big_integer::trivial_division::modulo(big_integer &dividend, const big_integer &divisor,
                                                   big_integer::multiplication_rule multiplication_rule) const
                                                   {
                                                        return dividend %= divisor;
                                                   }

big_integer &
big_integer::multiply(big_integer &first_multiplier, const big_integer &second_multiplier, allocator *allocator,
                      big_integer::multiplication_rule multiplication_rule)
{
   return first_multiplier *= second_multiplier;
}

big_integer
big_integer::multiply(const big_integer &first_multiplier, const big_integer &second_multiplier, allocator *allocator,
                      big_integer::multiplication_rule multiplication_rule)
{
    return big_integer(first_multiplier) *= second_multiplier;
}

big_integer &big_integer::divide(big_integer &dividend, const big_integer &divisor, allocator *allocator,
                                 big_integer::division_rule division_rule,
                                 big_integer::multiplication_rule multiplication_rule)
{
   return dividend /= divisor;
}

big_integer big_integer::divide(const big_integer &dividend, const big_integer &divisor, allocator *allocator,
                                big_integer::division_rule division_rule,
                                big_integer::multiplication_rule multiplication_rule)
{
   return big_integer(dividend) /= divisor;
}

big_integer &big_integer::modulo(big_integer &dividend, const big_integer &divisor, allocator *allocator,
                                 big_integer::division_rule division_rule,
                                 big_integer::multiplication_rule multiplication_rule)
{
  return  dividend %= divisor;
}

big_integer big_integer::modulo(const big_integer &dividend, const big_integer &divisor, allocator *allocator,
                                big_integer::division_rule division_rule,
                                big_integer::multiplication_rule multiplication_rule)
{
  return big_integer(dividend) %= divisor;
}