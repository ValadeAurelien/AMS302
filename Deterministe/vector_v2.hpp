#ifndef VECTOR_V2_H
#define VECTOR_V2_H
#include <vector>
#include <cmath>
#include <exception>
#include <utility>

/*
** Exception à renvoyer si les tailles ne 
** correspondent pas entre deux vectors
** sommés
*/
struct size_exception : std::exception
{
    const char* what() const noexcept { return "Size(s) do(es) not match !"; }
};

/*
** On créé ici une classe dérivée du vector sur 
** laquelle  on définit des opérations standards 
** telles que +, -, *, / et norm 
*/
template <typename T>
class vectorV2 : public std::vector<T>
{
public:
    vectorV2 () : std::vector<T> () {}
    vectorV2 (unsigned size) : std::vector<T> (size) {}
    vectorV2 (unsigned size, T val) : std::vector<T> (size, val) {}
    vectorV2 (const std::vector<T>& vec) : std::vector<T> (vec) {}

    vectorV2<T> operator+(const vectorV2& SV) const;
    vectorV2<T> operator-(const vectorV2& SV) const;
    vectorV2<T> operator*(T q) const;
    vectorV2<T> operator/(T q) const;
    double norm() const;
    void flip();
};

template <typename T>
vectorV2<T> vectorV2<T>::operator+(const vectorV2<T>& SV) const
{
    if (this->size() != SV.size()) throw size_exception();
    vectorV2<T> res(this->size());
    for (int i=0; i<this->size(); i++)
	res.at(i) = this->at(i) + SV.at(i);
    return res;
}

template <typename T>
vectorV2<T> vectorV2<T>::operator-(const vectorV2<T>& SV) const
{
    if (this->size() != SV.size()) throw size_exception();
    vectorV2<T> res(this->size());
    for (int i=0; i<this->size(); i++)
	res.at(i) = this->at(i) - SV.at(i);
    return res;
}

template<typename T>
vectorV2<T> vectorV2<T>::operator*(T q) const
{
    vectorV2<T> res(this->size());
    for (int i=0; i<this->size(); i++)
	res.at(i) = this->at(i) * q;
    return res;
}

template<typename T>
vectorV2<T> vectorV2<T>::operator/(T q) const
{
    vectorV2<T> res(this->size());
    for (int i=0; i<this->size(); i++)
	res.at(i) = this->at(i) / q;
    return res;
}

template<typename T>
double vectorV2<T>::norm() const
{
    double sum = 0;
    for (auto const& v : (*this)) sum += pow(v, 2);
    return sum;
}

template<typename T>
void vectorV2<T>::flip()
{
    size_t s = this->size();
    for (int i=0; i<s/2; ++i) 
	std::swap(this->at(i), this->at(s-1-i));
}

#endif
