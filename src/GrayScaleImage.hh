#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"

#pragma once


#include "Types.hh"
#include "Debug.hh"

#include <string>
#include <vector>


//*******************************************************************************************************************
/*! Class for loading and storing png images
*
* Very simple wrapper around lodepng library.
* You also have to compile and link lodepng.cc
*/
//*******************************************************************************************************************
class GrayScaleImage
{
public:
    /// Loads a grayscale png image from the specified file
    GrayScaleImage(const std::string &pngFilename);
    /// creates a grayscale png image with the given data
    GrayScaleImage(const std::string &pngFilename, unsigned int tmpWidth, unsigned int tmpHeight);

    void save(const std::string &pngFilename);

    GrayScaleImage getResizedImage(int newWidth, int newHeight) const;

    int width()  const
    {
        return size_[0];
    }
    int height() const
    {
        return size_[1];
    }

    int size(int coord) const;

    /// Returns a value between 0 and 1
    /// 0 means black - 1 means white
    real operator()(int x, int y) const;
    inline unsigned char &operator()(int x, int y);

    /// Returns the gray value of the specified pixel (between 0 and 255)
    unsigned char &getElement(int x, int y);
    unsigned char   getElement(int x, int y) const;



protected:
    GrayScaleImage() {} // required in getResizedImage
    int size_[2];                      //< 0=width,  1=height
    std::vector<unsigned char> image_; //< raw pixels
};







//===================================================================================================================
//
//  Implementation of inline functions
//
//===================================================================================================================

inline unsigned char &GrayScaleImage::operator()(int x, int y)
{
    int yFlip = size_[1] - y - 1;

    return getElement(x, yFlip);
}

inline unsigned char   &GrayScaleImage::getElement(int x, int y)
{
    ASSERT(x >= 0  && y >= 0);
    ASSERT(x < size_[0]);
    ASSERT(y < size_[1]);
    return image_[ y * size_[0] + x ];
}

inline unsigned char GrayScaleImage::getElement(int x, int y) const
{
    return const_cast<GrayScaleImage *>(this)->getElement(x, y);
}


inline int GrayScaleImage::size(int coord) const
{
    ASSERT(coord < 2);
    return size_[coord];
}


#pragma GCC diagnostic pop
