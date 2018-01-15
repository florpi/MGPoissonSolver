#include <iostream>
#include<vector>
using namespace std;

template <typename T>
class vector2d {
public:
    vector2d(size_t d1=0, size_t d2=0, T const & t=T()) :
        d1(d1), d2(d2),  data(d1*d2, t)
    {}

    T & operator()(size_t i, size_t j) {
        return data[i*d1 + j];
    }

    T const & operator()(size_t i, size_t j) const {
        return data[i*d1 + j];
    }
	void print();

private:
	size_t d1,d2;
	vector<T> data;
};

template <class T>
void vector2d<T>::print()
{
	for (size_t r = 0; r < d1; ++r)
	{
		for (size_t c = 0; c < d2; ++c)
		{
			cout << data[ (r*d1) + c ] << " ";
		}
		cout << endl;
	}
}
