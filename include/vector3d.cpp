#include <iostream>
#include<vector>
using namespace std;
template <typename T>
class vector3d {
public:
    vector3d(size_t d1=0, size_t d2=0, size_t d3=0, T const & t=T()) :
        d1(d1), d2(d2), d3(d3), data(d1*d2*d3, t)
    {}

    T & operator()(size_t i, size_t j, size_t k) {
        return data[i*d2*d3 + j*d3 + k];
    }

    T const & operator()(size_t i, size_t j, size_t k) const {
        return data[i*d2*d3 + j*d3 + k];
    }

private:
	size_t d1,d2,d3;
	vector<T> data;
};

/*
int main(){

	vector3d <double>  vec(3,3,3);
	vec(0,0,0) = 1.5;
	cout << vec(0,0,0) << endl;
	return 0;
}
*/
