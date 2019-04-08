#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}


int main()
{
    std::cerr << "In main!" << endl;

    // 1024x1350 img dimensions
    int imgWidth = 1024;
    int imgHeight = 1350;

    vtkImageData *image = NewImage(imgWidth, imgHeight);
    unsigned char *buffer =
      (unsigned char *) image->GetScalarPointer(0,0,0);

    // 27 strips, 50 pixels each strip
    int numStrips = 27;
    int stripHeight = 50;

    //red green and blue for each pixel in the strip
    int r = 0;
    int g = 0;
    int b = 0;

    int stripPixels = 3 * imgWidth * stripHeight;

    for(int i = 0; i < numStrips; i++){
      int red = i / 9;
      r =
        red == 0 ? 0
        : red == 1 ? 128
        : 255;

      int green = (i / 3) % 3;
      g =
        green == 0 ? 0
        : green == 1 ? 128
        : 255;

      int blue = i % 3;
      b =
        blue == 0 ? 0
        : blue == 1 ? 128
        : 255;

      for(int j = 0; j < stripPixels; j = j + 3){
        int imageOffset = i * stripPixels + j;
        buffer[imageOffset] = r;
        buffer[imageOffset + 1] = g;
        buffer[imageOffset + 2] = b;
      }
    }
    WriteImage(image, "proj1A");
}
