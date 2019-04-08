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
  vtkPolyDataReader *rdr = vtkPolyDataReader::New();
  rdr->SetFileName("proj1c_geometry.vtk");
  cerr << "Reading" << endl;
  rdr->Update();
  if (rdr->GetOutput()->GetNumberOfCells() == 0)
  {
      cerr << "Unable to open file!!" << endl;
      exit(EXIT_FAILURE);
  }
  vtkPolyData *pd = rdr->GetOutput();
  int numTris = pd->GetNumberOfCells();
  vtkPoints *pts = pd->GetPoints();
  vtkCellArray *cells = pd->GetPolys();
  vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
  float *color_ptr = colors->GetPointer(0);
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
      tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
      tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
      tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
  }
  cerr << "Done reading" << endl;

  return tris;
}

void placeDownwardTrianglePixels(Triangle triangle, unsigned char *buffer){
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
  maxY = baseY1 < 1343 ? floor_441(baseY1) : 1343;

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
    edgeX2 = pixelX2 < 1785 ? pixelX2 : 1785;

    //place colors along line
    for(int j = edgeX1; j <= edgeX2; j++){
      int imageOffset = 3*(i * 1786 + j);
      buffer[imageOffset] = triangle.color[0];
      buffer[imageOffset + 1] = triangle.color[1];
      buffer[imageOffset + 2] = triangle.color[2];
    }
  }
}

void placeUpwardTrianglePixels(Triangle triangle, unsigned char *buffer){
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
  minY = baseY1 > 0 ? ceil_441(baseY1) : 0;
  maxY = tipY < 1343 ? floor_441(tipY) : 1343;

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
      pixelX1 = ceil_441(baseX1 + ((double)i - baseY1)/slope1);
    }else{
      pixelX1 = ceil_441(baseX1);
    }
    if(slope2 != 0){
      pixelX2 = floor_441(baseX2 + (((double)i) - baseY2)/slope2);
    }else{
      pixelX2 = floor_441(baseX2);
    }

    //check endpoint for out of bounds
    int edgeX1, edgeX2;
    edgeX1 = pixelX1 > 0 ? pixelX1 : 0;
    edgeX2 = pixelX2 < 1785 ? pixelX2 : 1785;

    //place colors along line
    for(int j = edgeX1; j <= edgeX2; j++){
      int imageOffset = 3*(i * 1786 + j);
      buffer[imageOffset] = triangle.color[0];
      buffer[imageOffset + 1] = triangle.color[1];
      buffer[imageOffset + 2] = triangle.color[2];
    }
  }
}

//values of a triangle needed to create a pixel with 2
struct splitValues{
  double maxY, maxX, minY, minX, splitY, splitX, slope;
};

//find the min,max,seperation, and slope values
struct splitValues getValues(Triangle triangle){
  double minY = triangle.Y[0], minX = triangle.X[0];
  for(int i = 1; i<=2; i++){
    if(minY > triangle.Y[i]){
      minY = triangle.Y[i];
      minX = triangle.X[i];
    }
  }

  double maxY = triangle.Y[0], maxX = triangle.X[0];
  for(int i = 1; i<=2; i++){
    if(maxY < triangle.Y[i]){
      maxY = triangle.Y[i];
      maxX = triangle.X[i];
    }
  }

  double slope = (maxY - minY)/(maxX - minX);

  double splitY, splitX;
  for(int i = 0; i<=2; i++){
    if((minY < triangle.Y[i]) && (maxY > triangle.Y[i])){
      splitY = triangle.Y[i];
      splitX = triangle.X[i];
    }
  }

  return {maxY, maxX, minY, minX, splitY, splitX, slope};
}

void placeTrianglePixels(Triangle triangle,unsigned char *buffer){
  //check if triangle is up down or both
  if(triangle.Y[0] == triangle.Y[1]){
    if(triangle.Y[0] > triangle.Y[2]){
      placeDownwardTrianglePixels(triangle, buffer);
    }
    else{
      placeUpwardTrianglePixels(triangle, buffer);
    }
  }
  else if(triangle.Y[0] == triangle.Y[2]){
    if(triangle.Y[0] > triangle.Y[1]){
      placeDownwardTrianglePixels(triangle, buffer);
    }
    else{
      placeUpwardTrianglePixels(triangle, buffer);
    }
  }
  else if(triangle.Y[2] == triangle.Y[1]){
    if(triangle.Y[0] < triangle.Y[1]){
      placeDownwardTrianglePixels(triangle, buffer);
    }
    else{
      placeUpwardTrianglePixels(triangle, buffer);
    }
  }
  else{
    //find the points of intersection between the two triangles
    splitValues triangleValues = getValues(triangle);

    Triangle upwardTriangle, downwardTriangle;

    upwardTriangle.color[0] = downwardTriangle.color[0] = triangle.color[0];
    upwardTriangle.color[1] = downwardTriangle.color[1] = triangle.color[1];
    upwardTriangle.color[2] = downwardTriangle.color[2] = triangle.color[2];

    upwardTriangle.Y[1] = downwardTriangle.Y[1] = upwardTriangle.Y[2] = downwardTriangle.Y[2] =  triangleValues.splitY;
    upwardTriangle.X[1] = downwardTriangle.X[1] = triangleValues.splitX;
    upwardTriangle.X[2] = downwardTriangle.X[2] = (triangleValues.splitY - (triangleValues.maxY + -1 * triangleValues.slope * triangleValues.maxX))/triangleValues.slope;

    upwardTriangle.Y[0] = triangleValues.maxY;
    upwardTriangle.X[0] = triangleValues.maxX;
    downwardTriangle.Y[0] = triangleValues.minY;
    downwardTriangle.X[0] = triangleValues.minX;

    placeDownwardTrianglePixels(downwardTriangle, buffer);
    placeUpwardTrianglePixels(upwardTriangle, buffer);
  }
}

int main()
{
  vtkImageData *image = NewImage(1786, 1344);
  unsigned char *buffer =
    (unsigned char *) image->GetScalarPointer(0,0,0);
  int npixels = 1786*1344;
  for (int i = 0 ; i < npixels*3 ; i++)
      buffer[i] = 0;

  std::vector<Triangle> triangles = GetTriangles();

  Screen screen;
  screen.buffer = buffer;
  screen.width = 1786;
  screen.height = 1344;

  //places all traingles one by one
  for(int j = 0; j < triangles.size(); j++){
   placeTrianglePixels(triangles[j],buffer);
  }
  WriteImage(image, "project1C");
}
