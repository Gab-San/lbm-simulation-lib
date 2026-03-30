#pragma once

namespace lbm {

namespace utils {

template<typename T, int dim = 2>
struct Point;

template<typename T>
struct Point<T, 2> {
    const T x;
    const T y;
    Point(T _x, T _y) : x(_x), y(_y) {};
};


template<typename T>
struct Point<T, 3> {
    const T x;
    const T y;
    const T z;
    Point(T _x, T _y, T _z) : 
	x(_x), y(_y), z(_z) 
    {};
};

}

}
