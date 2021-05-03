#ifndef ARRAY_HH
#define ARRAY_HH

#include "Debug.hh"
#include "Types.hh"
#include <stdlib.h>
#include <iostream>
#include <string.h>


//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
template<typename T>
class ArrayTemplate
{
public:

    // Constructors for 1D,2D and 3D
    ArrayTemplate();
    ArrayTemplate(int xSize);
    ArrayTemplate(int xSize, int ySize);
    ArrayTemplate(int xSize, int ySize, int zSize);
    // Depending on your implementation you might need the following:
    ~ArrayTemplate();
    ArrayTemplate(const ArrayTemplate &s);


    // Access Operators for 1D, 2D and 3D
    inline T &operator()(int i);
    inline T &operator()(int i , int j);
    inline T &operator()(int i, int j, int k);


    // Assignment Operator
    inline ArrayTemplate<T> &operator = (const ArrayTemplate<T> &s);

    // for const Arrays the following access operators are required
    inline const T &operator()(int i) const;
    inline const T &operator()(int i , int j) const;
    inline const T &operator()(int i, int j, int k) const;

    // get max element of the array
    T maxE();

    // initialize the whole array with a constant value
    void fill(T value);


    // return total size of the array
    int getSize() const;

    // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
    // other dimension values are not allowed
    int getSize(int dimension) const;


    // Print the whole array ( for debugging purposes )
    void print();

private:

    // dimension sizes
    int x, y, z;

    // actual array
    T *arr;

};

typedef ArrayTemplate<real> Array;
typedef ArrayTemplate<flag> FlagArray;


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 1D
template<typename T>
inline T &ArrayTemplate<T>::operator()(int i)
{
    ASSERT_MSG((i >= 0), "wrong input for i:" << i);
    ASSERT_MSG((i < x), "wrong input for i:" << i << ", with maxlength = " << x);

    return arr[i];
}

// Operator() 2D
template<typename T>
inline T &ArrayTemplate<T>::operator()(int i, int j)
{
    ASSERT_MSG((i >= 0), "wrong input for i:" << i);
    ASSERT_MSG((i < x), "wrong input for i:" << i << ", with maxlength = " << x);
    ASSERT_MSG((j >= 0), "wrong input for j:" << j);
    ASSERT_MSG((j < y), "wrong input for j:" << j << ", with maxlength = " << y);

    return arr[i + j * x];
}

// Operator() 3D
template<typename T>
inline T &ArrayTemplate<T>::operator()(int i, int j, int k)
{
    ASSERT_MSG((i >= 0), "wrong input for i:" << i);
    ASSERT_MSG((i < x), "wrong input for i:" << i << ", with maxlength = " << x);
    ASSERT_MSG((j >= 0), "wrong input for j:" << j);
    ASSERT_MSG((j < y), "wrong input for j:" << j << ", with maxlength = " << y);
    ASSERT_MSG((k >= 0), "wrong input for k:" << k);
    ASSERT_MSG((k < z), "wrong input for k:" << k << ", with maxlength = " << z);

    return arr[i + j * x + k * x * y];
}

// Operator() 1D
template<typename T>
inline const T &ArrayTemplate<T>::operator()(int i) const
{
    ASSERT_MSG((i >= 0), "wrong input for i:" << i);
    ASSERT_MSG((i < x), "wrong input for i:" << i << ", with maxlength = " << x);

    return arr[i];
}

// Operator() 2D
template<typename T>
inline const T &ArrayTemplate<T>::operator()(int i, int j) const
{
    ASSERT_MSG((i >= 0), "wrong input for i:" << i);
    ASSERT_MSG((i < x), "wrong input for i:" << i << ", with maxlength = " << x);
    ASSERT_MSG((j >= 0), "wrong input for j:" << j);
    ASSERT_MSG((j < y), "wrong input for j:" << j << ", with maxlength = " << y);

    return arr[i + j * x];
}

// Operator() 3D
template<typename T>
inline const T &ArrayTemplate<T>::operator()(int i, int j, int k) const
{
    ASSERT_MSG((i >= 0), "wrong input for i:" << i);
    ASSERT_MSG((i < x), "wrong input for i:" << i << ", with maxlength = " << x);
    ASSERT_MSG((j >= 0), "wrong input for j:" << j);
    ASSERT_MSG((j < y), "wrong input for j:" << j << ", with maxlength = " << y);
    ASSERT_MSG((k >= 0), "wrong input for k:" << k);
    ASSERT_MSG((k < z), "wrong input for k:" << k << ", with maxlength = " << z);

    return arr[i + j * x + k * x * y];
}

template<typename T>
inline ArrayTemplate<T> &ArrayTemplate<T>::operator = (const ArrayTemplate<T> &s)
{
    PROG("assign array");
    CHECK_MSG(s.arr != NULL, "assignment with NULL"); // check s.arr

    // set dimension sizes and get memory
    x = s.getSize(0);
    y = s.getSize(1);
    z = s.getSize(2);

    if (arr != NULL)
        delete[] arr;

    arr = new T[s.getSize()];
    ASSERT_MSG(arr != NULL, "malloc failed"); // check if malloc was successful

    memcpy(arr, s.arr, s.getSize()*sizeof(T));   // copy the array

    return *this;
}

#endif //ARRAY_HHH
