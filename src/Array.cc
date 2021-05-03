#include "Array.hh"


//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

template<typename T>
ArrayTemplate<T>::ArrayTemplate()
    : x(0), y(0), z(0), arr(NULL) {}

template<typename T>
ArrayTemplate<T>::ArrayTemplate(int xSize)
{
    PROG("construct 1-dim array");
    CHECK_MSG((xSize >= 0), "wrong input for xSize: " << xSize);
    // set dimension sizes and get memory
    x = xSize;
    y = 0;
    z = 0;
    arr = new T[x];
    ASSERT_MSG(arr != NULL, "malloc failed"); // check if malloc was successful
}

template<typename T>
ArrayTemplate<T>::ArrayTemplate(int xSize, int ySize)
{
    PROG("construct 2-dim array");
    CHECK_MSG((xSize >= 0), "wrong input for xSize: " << xSize);
    CHECK_MSG((ySize >= 0), "wrong input for ySize: " << ySize);
    // set dimension sizes and get memory
    x = xSize;
    y = ySize;
    z = 0;
    arr = new T[x * y];
    ASSERT_MSG(arr != NULL, "malloc failed"); // check if malloc was successful
}

template<typename T>
ArrayTemplate<T>::ArrayTemplate(int xSize, int ySize, int zSize)
{
    PROG("construct 3-dim array");
    CHECK_MSG((xSize >= 0), "wrong input for xSize: " << xSize);
    CHECK_MSG((ySize >= 0), "wrong input for ySize: " << ySize);
    CHECK_MSG((zSize >= 0), "wrong input for zSize: " << zSize);
    // set dimension sizes and get memory
    x = xSize;
    y = ySize;
    z = zSize;
    arr = new T[x * y * z];
    ASSERT_MSG(arr != NULL, "malloc failed"); // check if malloc was successful
}

template<typename T>
ArrayTemplate<T>::~ArrayTemplate()
{
    // free memory
    if (arr != NULL)
        delete[] arr;
}

template<typename T>
ArrayTemplate<T>::ArrayTemplate(const ArrayTemplate &s)
{
    PROG("copy array") ;
    CHECK_MSG(s.arr != NULL, "copy with NULL"); // check s.arr
    // set dimension sizes and get memory
    x = s.getSize(0);
    y = s.getSize(1);
    z = s.getSize(2);
    arr = new T[s.getSize()];
    ASSERT_MSG(arr != NULL, "malloc failed"); // check if malloc was successful
    std::copy(s.arr, s.arr + s.getSize(), arr); // copy the array
}


//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================

// get max element of arr (without border)
template<typename T>
T ArrayTemplate<T>::maxE()
{
    T maximum = 0.0;
    T minimum = 0.0;

    if (y == 0)      // 1 dimension
    {
        for (int i = 1; i < x - 1; i++)
        {
            if (i == 1)
            {
                maximum = arr[i];
                minimum = arr[i];
            }
            if (maximum <= arr[i])
                maximum = arr[i];
            if (minimum >= arr[i])
                minimum = arr[i];
        }
    }

    if (z == 0)      // 2 dimensions
    {
        for (int i = 1; i < x - 1; i++)
        {
            for (int j = 1; j < y - 1; j++)
            {
                if (i == 1 && j == 1)
                {
                    maximum = arr[i + j * x];
                    minimum = arr[i + j * x];
                }
                if (maximum <= arr[i + j * x])
                    maximum = arr[i + j * x];
                if (minimum >= arr[i + j * x])
                    minimum = arr[i + j * x];
            }
        }
    }

    if (z != 0)      // 3 dimensions
    {
        for (int i = 1; i < x - 1; i++)
        {
            for (int j = 1; j < y - 1; j++)
            {
                for (int k = 1; k < z - 1; k++)
                {
                    if (i == 1 && j == 1 && k == 1)
                    {
                        maximum = arr[i + j * x + k * y];
                        minimum = arr[i + j * x + k * y];
                    }
                    if (maximum <= arr[i + j * x + k * y])
                        maximum = arr[i + j * x + k * y];
                    if (minimum >= arr[i + j * x + k * y])
                        minimum = arr[i + j * x + k * y];
                }
            }
        }
    }

    maximum = std::max(maximum, static_cast <T>((-1) * minimum));

    return maximum;
}

//initialize the whole array with a constant value
template<typename T>
void ArrayTemplate<T>::fill(T value)
{
    PROG("fill array with " << value);
    if (x != 0)   // for valid arrays
    {
        if (y == 0)   // 1 dimension
            std::fill(arr, arr + x, value); // fill the array

        if (z == 0)   // 2 dimensions
            std::fill(arr, arr + x * y, value); // fill the array
        if (z != 0)   // 3 dimensions
            std::fill(arr, arr + x * y * z, value); // fill the array
    }
}


// Print the whole array (for debugging purposes)
template<typename T>
void ArrayTemplate<T>::print()
{
    // For 2D Arrays the positive x-coordinate goes to the right
    //                   positive y-coordinate goes upwards
    //      -> the line with highest y-value should be printed first
    if (x != 0)     // for valid arrays
    {
        if (y == 0)    // 1 dimension
        {
            for (int i = 0; i < x; i++)
                std::cout << arr[i] << " ";

            std::cout << "\n\n";

            return;

        }

        if (z == 0)     // 2 dimensions
        {
            // print like in the explanation above
            for (int j = y - 1; j >= 0; j--)
            {
                for (int i = 0; i < x; i++)
                    std::cout << arr[i + j * x] << " ";

                std::cout << "\n";
            }
            std::cout << "\n" << std::endl;
        }
        else     // 3 dimensions
        {
            // print one 2-dim array after another
            for (int k = 0; k < z; k++)
            {
                for (int j = y - 1; j >= 0; j--)
                {
                    for (int i = 0; i < x; i++)
                        std::cout << arr[i + j * x + k * x * y] << " ";

                    std::cout << "\n";
                }
                std::cout << "\n" << std::endl;
            }
        }
    }
    else
    {
        std::cout << "NULL \n" << std::endl; // say "there's no array"
    }
}

template<typename T>
int ArrayTemplate<T>::getSize(int dimension) const
{
    // check input
    CHECK_MSG((dimension == 0 || dimension == 1 || dimension == 2), "Wrong input for dimension: " << dimension);

    if (x != 0)   // for valid arrays
    {
        if (dimension == 0) // 1 dimension
            return x;

        if (dimension == 1) // 2 dimensions
            return y;

        if (dimension == 2) // 3 dimensions
            return z;

    }

    return 0;
}

//return total size of the array
template<typename T>
int ArrayTemplate<T>::getSize() const
{
    if (x != 0)   // for valid arrays
    {
        if (y == 0)   // 1 dimension
            return x;

        if (z == 0)   // 2 dimensions
            return x * y;

        return x * y * z; // 3 dimensions
    }

    return 0;
}


template class ArrayTemplate<real>;
template class ArrayTemplate<flag>;