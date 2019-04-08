#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <algorithm>
#include <stdio.h>

using std::cerr;
using std::endl;

double max(double x, double y)
{
   return (x > y) ? x : y;
}

double min(double x, double y)
{
   return (x < y) ? x : y;
}

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
      double  X[3];
      double  Y[3];
      double  Z[3];
      double  color[3][3];
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
  vtkPolyDataReader *rdr = vtkPolyDataReader::New();
  rdr->SetFileName("proj1d_geometry.vtk");
  cerr << "Reading" << endl;
  rdr->Update();
  cerr << "Done reading" << endl;
  if (rdr->GetOutput()->GetNumberOfCells() == 0)
  {
      cerr << "Unable to open file!!" << endl;
      exit(EXIT_FAILURE);
  }
  vtkPolyData *pd = rdr->GetOutput();
  int numTris = pd->GetNumberOfCells();
  vtkPoints *pts = pd->GetPoints();
  vtkCellArray *cells = pd->GetPolys();
  vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
  float *color_ptr = var->GetPointer(0);
  std::vector<Triangle> tris(numTris);
  vtkIdType npts;
  vtkIdType *ptIds;
  int idx;
  for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
  {
      if (npts != 3)
      {
          cerr << "Non-triangles!! ???" << endl;
          exit(EXIT_FAILURE);
      }
      tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
      tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
      tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
      tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
      tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
      tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
      tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
      tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
      tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
      // 1->2 interpolate between light blue, dark blue
      // 2->2.5 interpolate between dark blue, cyan
      // 2.5->3 interpolate between cyan, green
      // 3->3.5 interpolate between green, yellow
      // 3.5->4 interpolate between yellow, orange
      // 4->5 interpolate between orange, brick
      // 5->6 interpolate between brick, salmon
      double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
      double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
      unsigned char RGB[8][3] = { { 71, 71, 219 },
                                  { 0, 0, 91 },
                                  { 0, 255, 255 },
                                  { 0, 128, 0 },
                                  { 255, 255, 0 },
                                  { 255, 96, 0 },
                                  { 107, 0, 0 },
                                  { 224, 76, 76 }
                                };
      for (int j = 0 ; j < 3 ; j++)
      {
          float val = color_ptr[ptIds[j]];
          int r;
          for (r = 0 ; r < 7 ; r++)
          {
              if (mins[r] <= val && val < maxs[r])
                  break;
          }
          if (r == 7)
          {
              cerr << "Could not interpolate color for " << val << endl;
              exit(EXIT_FAILURE);
          }
          double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
          tris[idx].color[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
          tris[idx].color[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
          tris[idx].color[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
      }
  }

  return tris;
}

void placeDownwardTrianglePixels(Triangle triangle, unsigned char *buffer, double* zvalues, int width, int height){
  //find which vertices are which
  double baseY1, baseX1, baseZ1, baseZ2, baseY2, baseX2, tipY, tipX, tipZ;
  double baseRed1, baseGreen1, baseBlue1, baseRed2, baseGreen2, baseBlue2, tipRed, tipGreen, tipBlue;
  if(triangle.Y[0] == triangle.Y[1]){
    baseY1 = triangle.Y[0];
    baseX1 = triangle.X[0] < triangle.X[1] ? triangle.X[0] : triangle.X[1];
    baseY2 = triangle.Y[0];
    baseX2 = triangle.X[0] < triangle.X[1] ? triangle.X[1] : triangle.X[0];
    baseZ1 = triangle.X[0] < triangle.X[1] ? triangle.Z[0]:triangle.Z[1];
    baseZ2 = triangle.X[0] < triangle.X[1] ? triangle.Z[1]:triangle.Z[0];
    tipY = triangle.Y[2];
    tipX = triangle.X[2];
    tipZ = triangle.Z[2];
    baseRed1 = triangle.X[0] < triangle.X[1] ? triangle.color[0][0] : triangle.color[1][0];
    baseGreen1 = triangle.X[0] < triangle.X[1] ? triangle.color[0][1] : triangle.color[1][1];
    baseBlue1 = triangle.X[0] < triangle.X[1] ? triangle.color[0][2] : triangle.color[1][2];
    baseRed2 = triangle.X[0] < triangle.X[1] ? triangle.color[1][0] : triangle.color[0][0];
    baseGreen2 = triangle.X[0] < triangle.X[1] ? triangle.color[1][1] : triangle.color[0][1];
    baseBlue2 = triangle.X[0] < triangle.X[1] ? triangle.color[1][2] : triangle.color[0][2];
    tipRed   = triangle.color[2][0];
    tipGreen = triangle.color[2][1];
    tipBlue  = triangle.color[2][2];
  }else if(triangle.Y[0] == triangle.Y[2]){
    baseY1 = triangle.Y[0];
    baseX1 = triangle.X[0] < triangle.X[2] ? triangle.X[0] : triangle.X[2];
    baseY2 = triangle.Y[0];
    baseX2 = triangle.X[0] < triangle.X[2] ? triangle.X[2] : triangle.X[0];
    baseZ1 = triangle.X[0] < triangle.X[2] ? triangle.Z[0]:triangle.Z[2];
    baseZ2 = triangle.X[0] < triangle.X[2] ? triangle.Z[2]:triangle.Z[0];
    tipY = triangle.Y[1];
    tipX = triangle.X[1];
    tipZ = triangle.Z[1];

    baseRed1 = triangle.X[0] < triangle.X[2] ? triangle.color[0][0] : triangle.color[2][0];
    baseGreen1 = triangle.X[0] < triangle.X[2] ? triangle.color[0][1] : triangle.color[2][1];
    baseBlue1 = triangle.X[0] < triangle.X[2] ? triangle.color[0][2] : triangle.color[2][2];
    baseRed2 = triangle.X[0] < triangle.X[2] ? triangle.color[2][0] : triangle.color[0][0];
    baseGreen2 = triangle.X[0] < triangle.X[2] ? triangle.color[2][1] : triangle.color[0][1];
    baseBlue2 = triangle.X[0] < triangle.X[2] ? triangle.color[2][2] : triangle.color[0][2];
    tipRed   = triangle.color[1][0];
    tipGreen = triangle.color[1][1];
    tipBlue  = triangle.color[1][2];
  }else{
    baseY1 = triangle.Y[1];
    baseX1 = triangle.X[1] < triangle.X[2] ? triangle.X[1] : triangle.X[2];
    baseY2 = triangle.Y[1];
    baseX2 = triangle.X[1] < triangle.X[2] ? triangle.X[2] : triangle.X[1];
    baseZ1 = triangle.X[1] < triangle.X[2] ? triangle.Z[1]:triangle.Z[2];
    baseZ2 = triangle.X[1] < triangle.X[2] ? triangle.Z[2]:triangle.Z[1];
    tipY = triangle.Y[0];
    tipX = triangle.X[0];
    tipZ = triangle.Z[0];

    baseRed1 = triangle.X[1] < triangle.X[2] ? triangle.color[1][0] : triangle.color[2][0];
    baseGreen1 = triangle.X[1] < triangle.X[2] ? triangle.color[1][1] : triangle.color[2][1];
    baseBlue1 = triangle.X[1] < triangle.X[2] ? triangle.color[1][2] : triangle.color[2][2];
    baseRed2 = triangle.X[1] < triangle.X[2] ? triangle.color[2][0] : triangle.color[1][0];
    baseGreen2 = triangle.X[1] < triangle.X[2] ? triangle.color[2][1] : triangle.color[1][1];
    baseBlue2 = triangle.X[1] < triangle.X[2] ? triangle.color[2][2] : triangle.color[1][2];
    tipRed   = triangle.color[0][0];
    tipGreen = triangle.color[0][1];
    tipBlue  = triangle.color[0][2];
  }

  //check if vertices are out of bounds
  int minY, maxY;
  minY = tipY >= 0 ? ceil_441(tipY) : 0;
  maxY = baseY1 < height ? floor_441(baseY1) : height-1;

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
  double red1, green1, blue1, red2, green2, blue2;
  //scanline algorithm
  for(int i = minY; i <= maxY; i++){
    double pixelZ1, pixelZ2, pixelX1, pixelX2;
    pixelZ1 = tipZ != baseZ1 ? baseZ1 + ((tipZ - baseZ1)/(tipY - baseY1))*(i - baseY2) :  baseZ1;
    pixelZ2 = tipZ != baseZ2 ? baseZ2 + ((tipZ - baseZ2)/(tipY - baseY2))*(i - baseY2) : baseZ2;
    red1 = tipRed != baseRed1 ? baseRed1 + ((tipRed - baseRed1)/(tipY - baseY1))*(i - baseY2) : baseRed1;
    red2 = tipRed != baseRed2 ? baseRed2 + ((tipRed - baseRed2)/(tipY - baseY2))*(i - baseY2) : baseRed2;
    green1 = tipGreen != baseGreen1 ? baseGreen1 + ((tipGreen - baseGreen1)/(tipY - baseY1))*(i - baseY2) : baseGreen1;
    green2 = tipGreen != baseGreen2 ? baseGreen2 + ((tipGreen - baseGreen2)/(tipY - baseY2))*(i - baseY2) : baseGreen2;
    blue1 = tipBlue != baseBlue1 ? baseBlue1 + ((tipBlue - baseBlue1)/(tipY - baseY1))*(i - baseY2) : baseBlue1;
    blue2  = tipBlue != baseBlue2 ? baseBlue2 + ((tipBlue - baseBlue2)/(tipY - baseY2))*(i - baseY2) : baseBlue2;

    if(slope1 != 0){
      pixelX1 = tipX + (i - tipY)/slope1;
    }else{
      pixelX1 = baseX1;
    }
    if(slope2 != 0){
      pixelX2 = tipX + ((i) - tipY)/slope2;
    }else{
      pixelX2 = baseX2;
    }

    //check endpoint for out of bounds
    int edgeX1, edgeX2;
    edgeX1 = pixelX1 >= 0 ? ceil_441(pixelX1) : 0;
    edgeX2 = pixelX2 < width ? floor_441(pixelX2) : width-1;

    double pixelZ;
    //place colors along line
    for(int j = edgeX1; j <= edgeX2; j++){
      pixelZ = (j - pixelX1) * ((pixelZ2 - pixelZ1)/(pixelX2 - pixelX1)) + pixelZ1;
      int imageOffset = 3*(i * width + j);
      if((pixelZ <= 0.0) && (pixelZ >= -1.0) && (pixelZ >= zvalues[i * width + j])){
        zvalues[i * width + j] = pixelZ;
        buffer[imageOffset] = ceil_441(255*((j - pixelX1)*((red2 - red1)/(pixelX2 - pixelX1)) + red1));
        buffer[imageOffset + 1] = ceil_441(255*((j - pixelX1)*((green2 - green1)/(pixelX2 - pixelX1)) + green1));
        buffer[imageOffset + 2] = ceil_441(255*((j - pixelX1)*((blue2 - blue1)/(pixelX2 - pixelX1)) + blue1));
    }
  }
  }
}

void placeUpwardTrianglePixels(Triangle triangle, unsigned char *buffer, double* zvalues, int width, int height){
  double baseY1, baseX1, baseZ1, baseZ2, baseY2, baseX2, tipY, tipX, tipZ;
  double baseRed1, baseGreen1, baseBlue1, baseRed2, baseGreen2, baseBlue2, tipRed, tipGreen, tipBlue;
  if(triangle.Y[0] == triangle.Y[1]){
    baseY1 = triangle.Y[0];
    baseX1 = triangle.X[0] < triangle.X[1] ? triangle.X[0] : triangle.X[1];
    baseY2 = triangle.Y[0];
    baseX2 = triangle.X[0] < triangle.X[1] ? triangle.X[1] : triangle.X[0];
    baseZ1 = triangle.X[0] < triangle.X[1] ? triangle.Z[0]:triangle.Z[1];
    baseZ2 = triangle.X[0] < triangle.X[1] ? triangle.Z[1]:triangle.Z[0];
    tipY = triangle.Y[2];
    tipX = triangle.X[2];
    tipZ = triangle.Z[2];
    baseRed1 = triangle.X[0] < triangle.X[1] ? triangle.color[0][0] : triangle.color[1][0];
    baseGreen1 = triangle.X[0] < triangle.X[1] ? triangle.color[0][1] : triangle.color[1][1];
    baseBlue1 = triangle.X[0] < triangle.X[1] ? triangle.color[0][2] : triangle.color[1][2];
    baseRed2 = triangle.X[0] < triangle.X[1] ? triangle.color[1][0] : triangle.color[0][0];
    baseGreen2 = triangle.X[0] < triangle.X[1] ? triangle.color[1][1] : triangle.color[0][1];
    baseBlue2 = triangle.X[0] < triangle.X[1] ? triangle.color[1][2] : triangle.color[0][2];
    tipRed   = triangle.color[2][0];
    tipGreen = triangle.color[2][1];
    tipBlue  = triangle.color[2][2];
  }else if(triangle.Y[0] == triangle.Y[2]){
    baseY1 = triangle.Y[0];
    baseX1 = triangle.X[0] < triangle.X[2] ? triangle.X[0] : triangle.X[2];
    baseY2 = triangle.Y[0];
    baseX2 = triangle.X[0] < triangle.X[2] ? triangle.X[2] : triangle.X[0];
    baseZ1 = triangle.X[0] < triangle.X[2] ? triangle.Z[0]:triangle.Z[2];
    baseZ2 = triangle.X[0] < triangle.X[2] ? triangle.Z[2]:triangle.Z[0];
    tipY = triangle.Y[1];
    tipX = triangle.X[1];
    tipZ = triangle.Z[1];

    baseRed1 = triangle.X[0] < triangle.X[2] ? triangle.color[0][0] : triangle.color[2][0];
    baseGreen1 = triangle.X[0] < triangle.X[2] ? triangle.color[0][1] : triangle.color[2][1];
    baseBlue1 = triangle.X[0] < triangle.X[2] ? triangle.color[0][2] : triangle.color[2][2];
    baseRed2 = triangle.X[0] < triangle.X[2] ? triangle.color[2][0] : triangle.color[0][0];
    baseGreen2 = triangle.X[0] < triangle.X[2] ? triangle.color[2][1] : triangle.color[0][1];
    baseBlue2 = triangle.X[0] < triangle.X[2] ? triangle.color[2][2] : triangle.color[0][2];
    tipRed   = triangle.color[1][0];
    tipGreen = triangle.color[1][1];
    tipBlue  = triangle.color[1][2];
  }else{
    baseY1 = triangle.Y[1];
    baseX1 = triangle.X[1] < triangle.X[2] ? triangle.X[1] : triangle.X[2];
    baseY2 = triangle.Y[1];
    baseX2 = triangle.X[1] < triangle.X[2] ? triangle.X[2] : triangle.X[1];
    baseZ1 = triangle.X[1] < triangle.X[2] ? triangle.Z[1]:triangle.Z[2];
    baseZ2 = triangle.X[1] < triangle.X[2] ? triangle.Z[2]:triangle.Z[1];
    tipY = triangle.Y[0];
    tipX = triangle.X[0];
    tipZ = triangle.Z[0];

    baseRed1 = triangle.X[1] < triangle.X[2] ? triangle.color[1][0] : triangle.color[2][0];
    baseGreen1 = triangle.X[1] < triangle.X[2] ? triangle.color[1][1] : triangle.color[2][1];
    baseBlue1 = triangle.X[1] < triangle.X[2] ? triangle.color[1][2] : triangle.color[2][2];
    baseRed2 = triangle.X[1] < triangle.X[2] ? triangle.color[2][0] : triangle.color[1][0];
    baseGreen2 = triangle.X[1] < triangle.X[2] ? triangle.color[2][1] : triangle.color[1][1];
    baseBlue2 = triangle.X[1] < triangle.X[2] ? triangle.color[2][2] : triangle.color[1][2];
    tipRed   = triangle.color[0][0];
    tipGreen = triangle.color[0][1];
    tipBlue  = triangle.color[0][2];
  }

  //check if vertices are out of bounds
  int minY, maxY;
  minY = baseY1 >= 0 ? ceil_441(baseY1) : 0;
  maxY = tipY < height ? floor_441(tipY) : height-1;

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
  double red1, green1, blue1, red2, green2, blue2;
  //scanline algorithm
  for(int i = minY; i <= maxY; i++){
    double pixelZ1, pixelZ2, pixelX1, pixelX2;
    pixelZ1 = tipZ != baseZ1 ? baseZ1 + ((tipZ - baseZ1)/(tipY - baseY1))*(i - baseY2) :  baseZ1;
    pixelZ2 = tipZ != baseZ2 ? baseZ2 + ((tipZ - baseZ2)/(tipY - baseY2))*(i - baseY2) : baseZ2;
    red1 = tipRed != baseRed1 ? baseRed1 + ((tipRed - baseRed1)/(tipY - baseY1))*(i - baseY2) : baseRed1;
    red2 = tipRed != baseRed2 ? baseRed2 + ((tipRed - baseRed2)/(tipY - baseY2))*(i - baseY2) : baseRed2;
    green1 = tipGreen != baseGreen1 ? baseGreen1 + ((tipGreen - baseGreen1)/(tipY - baseY1))*(i - baseY2) : baseGreen1;
    green2 = tipGreen != baseGreen2 ? baseGreen2 + ((tipGreen - baseGreen2)/(tipY - baseY2))*(i - baseY2) : baseGreen2;
    blue1 = tipBlue != baseBlue1 ? baseBlue1 + ((tipBlue - baseBlue1)/(tipY - baseY1))*(i - baseY2) : baseBlue1;
    blue2  = tipBlue != baseBlue2 ? baseBlue2 + ((tipBlue - baseBlue2)/(tipY - baseY2))*(i - baseY2) : baseBlue2;
    //find lines endpoints
    if(slope1 != 0){
      pixelX1 = baseX1 + (i - baseY1)/slope1;
    }else{
      pixelX1 = (baseX1);
    }
    if(slope2 != 0){
      pixelX2 = baseX2 + ((i) - baseY2)/slope2;
    }else{
      pixelX2 = (baseX2);
    }

    //check endpoint for out of bounds
    int edgeX1, edgeX2;
    edgeX1 = pixelX1 >= 0 ? ceil_441(pixelX1) : 0;
    edgeX2 = pixelX2 < width ? floor_441(pixelX2) : width-1;

    double pixelZ;
    //place colors along line
    for(int j = edgeX1; j <= edgeX2; j++){
      pixelZ = (j - pixelX1) * ((pixelZ2 - pixelZ1)/(pixelX2 - pixelX1)) + pixelZ1;
      int imageOffset = 3*(i * width + j);
      if((pixelZ <= 0.0) && (pixelZ >= -1.0) && (pixelZ >= zvalues[i * width + j])){
        zvalues[i * width + j] = pixelZ;
        buffer[imageOffset] = ceil_441(255*((j - pixelX1)*((red2 - red1)/(pixelX2 - pixelX1)) + red1));
        buffer[imageOffset + 1] = ceil_441(255*((j - pixelX1)*((green2 - green1)/(pixelX2 - pixelX1)) + green1));
        buffer[imageOffset + 2] = ceil_441(255*((j - pixelX1)*((blue2 - blue1)/(pixelX2 - pixelX1)) + blue1));
    }
  }
}
}

void placeTrianglePixels(Triangle triangle, unsigned char *buffer, double* zvalues, int width, int height){
  //check if triangle is up down or both
  if(triangle.Y[0] == triangle.Y[1]){
    if(triangle.Y[0] > triangle.Y[2]){
      placeDownwardTrianglePixels(triangle, buffer, zvalues, width, height);
    }
    else{
      placeUpwardTrianglePixels(triangle, buffer, zvalues, width, height);
    }
  }
  else if(triangle.Y[0] == triangle.Y[2]){
    if(triangle.Y[0] > triangle.Y[1]){
      placeDownwardTrianglePixels(triangle, buffer, zvalues, width, height);
    }
    else{
      placeUpwardTrianglePixels(triangle, buffer, zvalues, width, height);
    }
  }
  else if(triangle.Y[2] == triangle.Y[1]){
    if(triangle.Y[0] < triangle.Y[1]){
      placeDownwardTrianglePixels(triangle, buffer, zvalues, width, height);
    }
    else{
      placeUpwardTrianglePixels(triangle, buffer, zvalues, width, height);
    }
  }
  else{
    double minY, minX, minZ,minR, minG, minB;
    double maxY, maxX,maxZ, maxR, maxG,maxB;
    double splitY, splitX, splitZ, splitR, splitG, splitB;

    minY = maxY = splitY = triangle.Y[0];
    minX = maxX = splitX = triangle.X[0];
    minZ = maxZ = splitZ = triangle.Z[0];
    maxR = minR = splitR = triangle.color[0][0];
    minG = maxG = splitG = triangle.color[0][1];
    minB = maxB = splitB = triangle.color[0][2];

    for(int i = 1; i<=2; i++){
      if(maxY < triangle.Y[i]){
        maxY = triangle.Y[i];
        maxX = triangle.X[i];
        maxZ = triangle.Z[i];
        maxR = triangle.color[i][0];
        maxG = triangle.color[i][1];
        maxB = triangle.color[i][2];
      }
    }
    for(int i = 1; i<=2; i++){
      if(minY > triangle.Y[i]){
        minY = triangle.Y[i];
        minX = triangle.X[i];
        minZ = triangle.Z[i];
        minR = triangle.color[i][0];
        minG = triangle.color[i][1];
        minB = triangle.color[i][2];
      }
    }
    for(int i = 0; i<=2; i++){
      if((maxY > triangle.Y[i]) && (minY < triangle.Y[i])){
        splitY = triangle.Y[i];
        splitX = triangle.X[i];
        splitZ = triangle.Z[i];
        splitR = triangle.color[i][0];
        splitG = triangle.color[i][1];
        splitB = triangle.color[i][2];
      }
    }

    Triangle upwardTriangle,downwardTriangle;

    upwardTriangle.Y[0] = maxY;
    upwardTriangle.X[0] = maxX;
    upwardTriangle.Z[0] = maxZ;
    downwardTriangle.Y[0] = minY;
    downwardTriangle.X[0] = minX;
    downwardTriangle.Z[0] = minZ;

    upwardTriangle.Y[1] = downwardTriangle.Y[1] = splitY;
    upwardTriangle.X[1] = downwardTriangle.X[1] = splitX;
    upwardTriangle.Z[1] = downwardTriangle.Z[1] = splitZ;

    upwardTriangle.Y[2] = downwardTriangle.Y[2] = splitY;
    upwardTriangle.X[2] = downwardTriangle.X[2] = minX + ((maxX - minX)/(maxY - minY))*(splitY - minY);
    upwardTriangle.Z[2] = downwardTriangle.Z[2] = minZ + ((maxZ - minZ)/(maxY - minY))*(splitY - minY);

    upwardTriangle.color[0][0] = maxR;
    upwardTriangle.color[0][1] = maxG;
    upwardTriangle.color[0][2] = maxB;

    downwardTriangle.color[0][0] = minR;
    downwardTriangle.color[0][1] = minG;
    downwardTriangle.color[0][2] = minB;

    upwardTriangle.color[1][0] = downwardTriangle.color[1][0] = splitR;
    upwardTriangle.color[1][1] = downwardTriangle.color[1][1] = splitG;
    upwardTriangle.color[1][2] = downwardTriangle.color[1][2] = splitB;

    upwardTriangle.color[2][0] = downwardTriangle.color[2][0] = minR + ((maxR - minR)/(maxY - minY))*(splitY - minY);
    upwardTriangle.color[2][1] = downwardTriangle.color[2][1] = minG + ((maxG - minG)/(maxY - minY))*(splitY - minY);
    upwardTriangle.color[2][2] = downwardTriangle.color[2][2] = minB + ((maxB - minB)/(maxY - minY))*(splitY - minY);

    placeDownwardTrianglePixels(downwardTriangle, buffer, zvalues, width, height);
    placeUpwardTrianglePixels(upwardTriangle, buffer, zvalues, width, height);

  }
}

int main()
{
  int width = 1000, height = 1000;
  vtkImageData *image = NewImage(width, height);
  unsigned char *buffer =
    (unsigned char *) image->GetScalarPointer(0,0,0);
  int npixels = width*height;
  double zvalues[npixels];
  std::fill_n(zvalues, npixels, -1.0);
  for (int i = 0 ; i < npixels*3 ; i++)
      buffer[i] = 0;

  std::vector<Triangle> triangles = GetTriangles();

  Screen screen;
  screen.buffer = buffer;
  screen.width = width;
  screen.height = height;

  //places all traingles one by one
  for(int j = 0; j < triangles.size(); j++){
   placeTrianglePixels(triangles[j],buffer,zvalues,width,height);
  }
  WriteImage(image, "allTriangles");
}
