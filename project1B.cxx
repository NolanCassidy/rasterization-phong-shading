#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}


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

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      unsigned char color[3];

  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;

  // would some methods for accessing and setting pixels be helpful?
};

std::vector<Triangle>
GetTriangles(void)
{
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204},
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       if (i == 50)
           rv[i].X[firstPt] = -10;
       rv[i].Y[firstPt] = posJ+10*(idxJ+1);
       rv[i].X[(firstPt+1)%3] = posI+99;
       rv[i].Y[(firstPt+1)%3] = posJ+10*(idxJ+1);
       rv[i].X[(firstPt+2)%3] = posI+i;
       rv[i].Y[(firstPt+2)%3] = posJ;
       if (i == 5)
          rv[i].Y[(firstPt+2)%3] = -50;
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;
}

void placeTrianglePixels(Triangle triangle, unsigned char *buffer){
  //determine which vertices are which
  double baseY1, baseX1, baseY2, baseX2, tipY, tipX;
  if(triangle.Y[0] == triangle.Y[1]){
    baseY1 = triangle.Y[0];
    baseX1 = triangle.X[0] < triangle.X[1] ? triangle.X[0] : triangle.X[1];
    baseY2 = triangle.Y[0];
    baseX2 = triangle.X[0] < triangle.X[1] ? triangle.X[1] : triangle.X[0];
    tipY = triangle.Y[2];
    tipX = triangle.X[2];
  }else if(triangle.Y[0] == triangle.Y[2]){
    baseY1 = triangle.Y[0];
    baseX1 = triangle.X[0] < triangle.X[2] ? triangle.X[0] : triangle.X[2];
    baseY2 = triangle.Y[0];
    baseX2 = triangle.X[0] < triangle.X[2] ? triangle.X[2] : triangle.X[0];
    tipY = triangle.Y[1];
    tipX = triangle.X[1];
  }else{
    baseY1 = triangle.Y[1];
    baseX1 = triangle.X[1] < triangle.X[2] ? triangle.X[1] : triangle.X[2];
    baseY2 = triangle.Y[1];
    baseX2 = triangle.X[1] < triangle.X[2] ? triangle.X[2] : triangle.X[1];
    tipY = triangle.Y[0];
    tipX = triangle.X[0];
  }

  //check if vertices are out of bounds
  int minY, maxY;
  minY = tipY > 0 ? ceil_441(tipY) : 0;
  maxY = baseY1 < 999 ? floor_441(baseY1) : 999;

  //find edge slopes
  double slope1, slope2;
  if(baseX1 != tipX){
    slope1 = (tipY - baseY1)/(tipX - baseX1);
  }else{
    slope1 = 0;
  }
  if(baseX2 != tipX){
    slope2 = (baseY2  - tipY)/(baseX2  - tipX);
  }else{
    slope2 = 0;
  }

  //scanline algorithm
  for(int i = minY; i <= maxY; i++){
    //find lines endpoints
    int pixelX1, pixelX2;
    if(slope1 != 0){
      pixelX1 = ceil_441(tipX + ((double)i - tipY)/slope1);
    }else{
      pixelX1 = ceil_441(baseX1);
    }
    if(slope2 != 0){
      pixelX2 = floor_441(tipX + (((double)i) - tipY)/slope2);
    }else{
      pixelX2 = floor_441(baseX2);
    }

    //check endpoint for out of bounds
    int edgeX1, edgeX2;
    edgeX1 = pixelX1 > 0 ? pixelX1 : 0;
    edgeX2 = pixelX2 < 999 ? pixelX2 : 999;

    //place colors along line
    for(int j = edgeX1; j <= edgeX2; j++){
      int imageOffset = 3*(i * 1000 + j);
      buffer[imageOffset] = triangle.color[0];
      buffer[imageOffset + 1] = triangle.color[1];
      buffer[imageOffset + 2] = triangle.color[2];
    }
  }
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);
   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;

    //place 100 triangles into buffer
    for(int j = 0; j<=99; j++){
      placeTrianglePixels(triangles[j], buffer);
    }
    WriteImage(image, "allTriangles");
}
