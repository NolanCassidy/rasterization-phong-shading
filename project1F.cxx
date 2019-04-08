#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#define NORMALS
using std::cerr;
using std::endl;

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };


    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;
double angleView[3];

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

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

double max(double x, double y, double z)
{
  return  (x < y) ? max(y,z) : max(x,z);
}

double min(double x, double y, double z)
{
  return  (x > y) ? min(y,z) : min(x,z);
}

vtkImageData * NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

	return img;
}

void WriteImage(vtkImageData *img, const char *filename){
    std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Screen{
    public:
        unsigned char   *buffer;
        int width, height;
};

double getDotProduct(double m1[], double m2[], int size){
    double dot = 0;
    for (int i=0; i < size; i++){ dot += (m1[i] * m2[i]); }
    return dot;
}

std::vector<double> getCrossProduct(double * m1, double * m2){
    std::vector<double> cross(3);
    cross[0] = (m1[1] * m2[2]) - (m1[2] * m2[1]);
    cross[1] = (m1[2] * m2[0]) - (m1[0] * m2[2]);
    cross[2] = (m1[0] * m2[1]) - (m1[1] * m2[0]);
    return cross;
}

double phongShade(LightingParameters lp, double *angleView, double *normal){
  double angle[3], R[3];
  for(int i = 0; i<=2; i++){
    angle[i] = 2*normal[i]*getDotProduct(lp.lightDir, normal, 3);
  }
  for(int i = 0; i<=2; i++){
    R[i] = angle[i] - lp.lightDir[i];
  }

  double diffuse,specular, shade;
  diffuse  = fabs(getDotProduct(lp.lightDir, normal, 3));
  specular = fmax(0, pow(getDotProduct(R, angleView,3), lp.alpha));
  shade = lp.Ka + lp.Kd*diffuse + lp.Ks*specular;

  return shade;
}

class Camera{
  public:
        double near, far;
        double angle;
        double position[3];
        double focus[3];
        double up[3];


        Matrix ViewTransform;
   	    Matrix CameraTransform;
        Matrix DeviceTransform;

        double O[3];
        double v1[3];
        double v2[3];
        double v3[3];

        void viewTransform(){
                  for(int i = 0; i<=3; i++){
                      for(int j = 0; j<=3; j++){
                        ViewTransform.A[i][j] = 0;
                      }
                  }

                  ViewTransform.A[0][0] = ViewTransform.A[1][1] = 1 / tan(angle/2);
                  ViewTransform.A[2][2] = (far + near) / (far - near);
                  ViewTransform.A[2][3] = -1;
                  ViewTransform.A[3][2] = (2 * far * near) / (far - near);

              }

        void cameraTransformHelper(double v[3]){
          double count = 0.0;
          for (int i = 0; i <=2; i++) {
              count += v[i] * v[i];
          }
          for (int i = 0; i <=2; i++) {
              v[i] = fabs(count) > 0.00001? v[i] / sqrt(count) : 0;
          }
        }

      	void cameraTransform(){
                  double f[3];
                  for (int i = 0; i <= 2; i++){ O[i] = position[i]; }
                  for (int i = 0; i <= 2; i++){ f[i] = O[i] - focus[i]; }

                  for (int i = 0; i <= 2; i++){ v1[i] = getCrossProduct(up, f)[i]; }
                  cameraTransformHelper(v1);

                  for (int i = 0; i <=2; i++){ v2[i] = getCrossProduct(f, v1)[i]; }
                  cameraTransformHelper(v2);

                  for (int i = 0; i <=2; i++){ v3[i] = f[i]; }
                  cameraTransformHelper(v3);

                  CameraTransform.A[0][3] = CameraTransform.A[1][3] = CameraTransform.A[2][3] = 0;
                  CameraTransform.A[3][3] = 1;

                  double iO[3];
                  for (int i = 0; i <=2; i++){ iO[i] = 0 - O[i]; }

                  CameraTransform.A[3][0] = getDotProduct(v1, iO, 3);
                  CameraTransform.A[3][1] = getDotProduct(v2, iO, 3);
                  CameraTransform.A[3][2] = getDotProduct(v3, iO, 3);

                  CameraTransform.A[0][0] = v1[0];
                  CameraTransform.A[1][0] = v1[1];
                  CameraTransform.A[2][0] = v1[2];

                  CameraTransform.A[0][1] = v2[0];
                  CameraTransform.A[1][1] = v2[1];
                  CameraTransform.A[2][1] = v2[2];

                  CameraTransform.A[0][2] = v3[0];
                  CameraTransform.A[1][2] = v3[1];
                  CameraTransform.A[2][2] = v3[2];
              }

        void deviceTransform(Screen screen){
            for(int i = 0; i<=3; i++){
                for(int j = 0; j<=3; j++){
                  DeviceTransform.A[i][j] = 0;
                }
            }
            DeviceTransform.A[0][0] = DeviceTransform.A[1][1] = DeviceTransform.A[3][0] = DeviceTransform.A[3][1] = screen.width / 2;
            DeviceTransform.A[2][2] = DeviceTransform.A[3][3] = 1;
        }
};

double SineParameterize(int curFrame, int nFrames, int ramp){
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera GetCamera(int frame, int nframes){
    double iO = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*iO);
    c.position[1] = 40*cos(2*M_PI*iO);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

class Triangle{
    public:
        double X[3];
        double Y[3];
        double Z[3];
        double color[3][3];
        double normals[3][3];
};

std::vector<Triangle> GetTriangles(void){
  vtkPolyDataReader *rdr = vtkPolyDataReader::New();
  rdr->SetFileName("proj1e_geometry.vtk");
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
  vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
  double *color_ptr = var->GetPointer(0);
  //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
  //float *color_ptr = var->GetPointer(0);
  vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
  float *normals = n->GetPointer(0);
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
      double *pt = NULL;
      pt = pts->GetPoint(ptIds[0]);
      tris[idx].X[0] = pt[0];
      tris[idx].Y[0] = pt[1];
      tris[idx].Z[0] = pt[2];
#ifdef NORMALS
      tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
      tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
      tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
      pt = pts->GetPoint(ptIds[1]);
      tris[idx].X[1] = pt[0];
      tris[idx].Y[1] = pt[1];
      tris[idx].Z[1] = pt[2];
#ifdef NORMALS
      tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
      tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
      tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
      pt = pts->GetPoint(ptIds[2]);
      tris[idx].X[2] = pt[0];
      tris[idx].Y[2] = pt[1];
      tris[idx].Z[2] = pt[2];
#ifdef NORMALS
      tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
      tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
      tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

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

Matrix getObjectMatrix(Camera camera,Screen s){
  //s1->s2->s3
  return Matrix::ComposeMatrices(Matrix::ComposeMatrices(camera.CameraTransform, camera.ViewTransform), camera.DeviceTransform);
}

Triangle TriangleTransform(Triangle triangle, Matrix matrix){
    Triangle t;

    //get new points for trianlge
    for (int i = 0; i <= 2; i++)
    {
        double oldT[4];
        double newT[4];
        oldT[0] = triangle.X[i];
        oldT[1] = triangle.Y[i];
        oldT[2] = triangle.Z[i];
        oldT[3] = 1;

        matrix.TransformPoint(oldT, newT);

        t.X[i] = newT[0] / newT[3];
        t.Y[i] = newT[1] / newT[3];
        t.Z[i] = newT[2] / newT[3];
    }

    //copy color over
    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2; j++) {
            t.color[i][j] = triangle.color[i][j];
            t.normals[i][j] = triangle.normals[i][j];
        }
    }

    return t;
}

void placeDownwardTrianglePixels(Triangle triangle, unsigned char *buffer, double* zvalues, int width, int height){
  //find which vertices are which
  double baseY1, baseX1, baseZ1, baseZ2, baseY2, baseX2, tipY, tipX, tipZ;
  double baseRed1, baseGreen1, baseBlue1, baseRed2, baseGreen2, baseBlue2, tipRed, tipGreen, tipBlue;
  double baseXNorm1, baseYNorm1, baseZNorm1, baseXNorm2, baseYNorm2, baseZNorm2, tipXNorm, tipYNorm, tipZNorm;
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
    baseXNorm1 = triangle.X[0] < triangle.X[1] ? triangle.normals[0][0] : triangle.normals[1][0];
    baseYNorm1 = triangle.X[0] < triangle.X[1] ? triangle.normals[0][1] : triangle.normals[1][1];
    baseZNorm1 = triangle.X[0] < triangle.X[1] ? triangle.normals[0][2] : triangle.normals[1][2];
    baseXNorm2 = triangle.X[0] < triangle.X[1] ? triangle.normals[1][0] : triangle.normals[0][0];
    baseYNorm2 = triangle.X[0] < triangle.X[1] ? triangle.normals[1][1] : triangle.normals[0][1];
    baseZNorm2 = triangle.X[0] < triangle.X[1] ? triangle.normals[1][2] : triangle.normals[0][2];
    tipXNorm   = triangle.normals[2][0];
    tipYNorm = triangle.normals[2][1];
    tipZNorm  = triangle.normals[2][2];
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
    baseXNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[0][0] : triangle.normals[2][0];
    baseYNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[0][1] : triangle.normals[2][1];
    baseZNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[0][2] : triangle.normals[2][2];
    baseXNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][0] : triangle.normals[0][0];
    baseYNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][1] : triangle.normals[0][1];
    baseZNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][2] : triangle.normals[0][2];
    tipXNorm   = triangle.normals[1][0];
    tipYNorm = triangle.normals[1][1];
    tipZNorm  = triangle.normals[1][2];
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
    baseXNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[1][0] : triangle.normals[2][0];
    baseYNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[1][1] : triangle.normals[2][1];
    baseZNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[1][2] : triangle.normals[2][2];
    baseXNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][0] : triangle.normals[1][0];
    baseYNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][1] : triangle.normals[1][1];
    baseZNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][2] : triangle.normals[1][2];
    tipXNorm   = triangle.normals[0][0];
    tipYNorm = triangle.normals[0][1];
    tipZNorm  = triangle.normals[0][2];
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
  double normX1, normY1, normZ1, normX2, normY2, normZ2;

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
    normX1 = tipXNorm != baseXNorm1 ? baseXNorm1 + ((tipXNorm - baseXNorm1)/(tipY - baseY1))*(i - baseY2) : baseXNorm1;
    normX2 = tipXNorm != baseXNorm2 ? baseXNorm2 + ((tipXNorm - baseXNorm2)/(tipY - baseY2))*(i - baseY2) : baseXNorm2;
    normY1 = tipYNorm != baseYNorm1 ? baseYNorm1 + ((tipYNorm - baseYNorm1)/(tipY - baseY1))*(i - baseY2) : baseYNorm1;
    normY2 = tipYNorm != baseYNorm2 ? baseYNorm2 + ((tipYNorm - baseYNorm2)/(tipY - baseY2))*(i - baseY2) : baseYNorm2;
    normZ1 = tipZNorm != baseZNorm1 ? baseZNorm1 + ((tipZNorm - baseZNorm1)/(tipY - baseY1))*(i - baseY2) : baseZNorm1;
    normZ2  = tipZNorm != baseZNorm2 ? baseZNorm2 + ((tipZNorm - baseZNorm2)/(tipY - baseY2))*(i - baseY2) : baseZNorm2;

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
    //place color along line
    for(int j = edgeX1; j <= edgeX2; j++){
      double normal[3];
      normal[0] = ((normX2 - normX1)*(j - pixelX1)/(pixelX2 - pixelX1)) + normX1;
      normal[1] = ((normY2 - normY1)*(j - pixelX1)/(pixelX2 - pixelX1)) + normY1;
      normal[2] = ((normZ2 - normZ1)*(j - pixelX1)/(pixelX2 - pixelX1)) + normZ1;

      double shade = phongShade(lp, angleView, normal);

      pixelZ = ((double)j - pixelX1) * ((pixelZ2 - pixelZ1)/(pixelX2 - pixelX1)) + pixelZ1;
      int imageOffset = 3*(i * width + j);
      if((pixelZ <= 0.0) && (pixelZ >= -1.0) && (pixelZ >= zvalues[i * width + j])){
        zvalues[i * width + j] = pixelZ;
        buffer[imageOffset] = min(255,ceil_441(255*((j - pixelX1)*((red2 - red1)/(pixelX2 - pixelX1)) + red1)*shade));
        buffer[imageOffset + 1] = min(255,ceil_441(255*((j - pixelX1)*((green2 - green1)/(pixelX2 - pixelX1)) + green1)*shade));
        buffer[imageOffset + 2] = min(255,ceil_441(255*((j - pixelX1)*((blue2 - blue1)/(pixelX2 - pixelX1)) + blue1)*shade));
    }
  }
  }
}

void placeUpwardTrianglePixels(Triangle triangle, unsigned char *buffer, double* zvalues, int width, int height){
  double baseY1, baseX1, baseZ1, baseZ2, baseY2, baseX2, tipY, tipX, tipZ;
  double baseRed1, baseGreen1, baseBlue1, baseRed2, baseGreen2, baseBlue2, tipRed, tipGreen, tipBlue;
  double baseXNorm1, baseYNorm1, baseZNorm1, baseXNorm2, baseYNorm2, baseZNorm2, tipXNorm, tipYNorm, tipZNorm;
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
    baseXNorm1 = triangle.X[0] < triangle.X[1] ? triangle.normals[0][0] : triangle.normals[1][0];
    baseYNorm1 = triangle.X[0] < triangle.X[1] ? triangle.normals[0][1] : triangle.normals[1][1];
    baseZNorm1 = triangle.X[0] < triangle.X[1] ? triangle.normals[0][2] : triangle.normals[1][2];
    baseXNorm2 = triangle.X[0] < triangle.X[1] ? triangle.normals[1][0] : triangle.normals[0][0];
    baseYNorm2 = triangle.X[0] < triangle.X[1] ? triangle.normals[1][1] : triangle.normals[0][1];
    baseZNorm2 = triangle.X[0] < triangle.X[1] ? triangle.normals[1][2] : triangle.normals[0][2];
    tipXNorm   = triangle.normals[2][0];
    tipYNorm = triangle.normals[2][1];
    tipZNorm  = triangle.normals[2][2];
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
    baseXNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[0][0] : triangle.normals[2][0];
    baseYNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[0][1] : triangle.normals[2][1];
    baseZNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[0][2] : triangle.normals[2][2];
    baseXNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][0] : triangle.normals[0][0];
    baseYNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][1] : triangle.normals[0][1];
    baseZNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][2] : triangle.normals[0][2];
    tipXNorm   = triangle.normals[1][0];
    tipYNorm = triangle.normals[1][1];
    tipZNorm  = triangle.normals[1][2];
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
    baseXNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[1][0] : triangle.normals[2][0];
    baseYNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[1][1] : triangle.normals[2][1];
    baseZNorm1 = triangle.X[0] < triangle.X[2] ? triangle.normals[1][2] : triangle.normals[2][2];
    baseXNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][0] : triangle.normals[1][0];
    baseYNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][1] : triangle.normals[1][1];
    baseZNorm2 = triangle.X[0] < triangle.X[2] ? triangle.normals[2][2] : triangle.normals[1][2];
    tipXNorm   = triangle.normals[0][0];
    tipYNorm = triangle.normals[0][1];
    tipZNorm  = triangle.normals[0][2];
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
  double normX1, normX2, normY1, normY2, normZ1, normZ2;

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
    normX1 = tipXNorm != baseXNorm1 ? baseXNorm1 + ((tipXNorm - baseXNorm1)/(tipY - baseY1))*(i - baseY2) : baseXNorm1;
    normX2 = tipXNorm != baseXNorm2 ? baseXNorm2 + ((tipXNorm - baseXNorm2)/(tipY - baseY2))*(i - baseY2) : baseXNorm2;
    normY1 = tipYNorm != baseYNorm1 ? baseYNorm1 + ((tipYNorm - baseYNorm1)/(tipY - baseY1))*(i - baseY2) : baseYNorm1;
    normY2 = tipYNorm != baseYNorm2 ? baseYNorm2 + ((tipYNorm - baseYNorm2)/(tipY - baseY2))*(i - baseY2) : baseYNorm2;
    normZ1 = tipZNorm != baseZNorm1 ? baseZNorm1 + ((tipZNorm - baseZNorm1)/(tipY - baseY1))*(i - baseY2) : baseZNorm1;
    normZ2  = tipZNorm != baseZNorm2 ? baseZNorm2 + ((tipZNorm - baseZNorm2)/(tipY - baseY2))*(i - baseY2) : baseZNorm2;

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
    //place color along line
    for(int j = edgeX1; j <= edgeX2; j++){
      double normal[3];
      normal[0] = ((normX2 - normX1)*(j - pixelX1)/(pixelX2 - pixelX1)) + normX1;
      normal[1] = ((normY2 - normY1)*(j - pixelX1)/(pixelX2 - pixelX1)) + normY1;
      normal[2] = ((normZ2 - normZ1)*(j - pixelX1)/(pixelX2 - pixelX1)) + normZ1;

      double shade = phongShade(lp, angleView, normal);

      int imageOffset = 3*(i * width + j);
      pixelZ = (j - pixelX1) * ((pixelZ2 - pixelZ1)/(pixelX2 - pixelX1)) + pixelZ1;
      if((pixelZ <= 0.0) && (pixelZ >= -1.0) && (pixelZ >= zvalues[i * width + j])){
        zvalues[i * width + j] = pixelZ;
        buffer[imageOffset] = min(255,ceil_441(255*((j - pixelX1)*((red2 - red1)/(pixelX2 - pixelX1)) + red1)*shade));
        buffer[imageOffset + 1] = min(255,ceil_441(255*((j - pixelX1)*((green2 - green1)/(pixelX2 - pixelX1)) + green1)*shade));
        buffer[imageOffset + 2] = min(255,ceil_441(255*((j - pixelX1)*((blue2 - blue1)/(pixelX2 - pixelX1)) + blue1)*shade));
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
    double minNormX, minNormY, minNormZ;
    double maxNormX, maxNormY, maxNormZ;
    double splitNormX, splitNormY, splitNormZ;

    minY = maxY = splitY = triangle.Y[0];
    minX = maxX = splitX = triangle.X[0];
    minZ = maxZ = splitZ = triangle.Z[0];
    maxR = minR = splitR = triangle.color[0][0];
    minG = maxG = splitG = triangle.color[0][1];
    minB = maxB = splitB = triangle.color[0][2];
    minNormY = maxNormY = splitNormY = triangle.normals[0][1];
    minNormX = maxNormX = splitNormX = triangle.normals[0][0];
    minNormZ = maxNormZ = splitNormZ = triangle.normals[0][2];

    for(int i = 1; i<=2; i++){
      if(maxY < triangle.Y[i]){
        maxY = triangle.Y[i];
        maxX = triangle.X[i];
        maxZ = triangle.Z[i];
        maxR = triangle.color[i][0];
        maxG = triangle.color[i][1];
        maxB = triangle.color[i][2];
        maxNormY = triangle.normals[i][1];
        maxNormX = triangle.normals[i][0];
        maxNormZ = triangle.normals[i][2];
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
        minNormY = triangle.normals[i][1];
        minNormX = triangle.normals[i][0];
        minNormZ = triangle.normals[i][2];
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
        splitNormY = triangle.normals[i][1];
        splitNormX = triangle.normals[i][0];
        splitNormZ = triangle.normals[i][2];
      }
    }
    double newX, newZ, newRed, newGreen, newBlue, newNormX, newNormY, newNormZ;

    //determine the x, z, and rgb values for the new point, generated by splitting the triangle in two, through interpolation
    newX = minX + ((maxX - minX)/(maxY - minY))*(splitY - minY);
    newZ = minZ + ((maxZ - minZ)/(maxY - minY))*(splitY - minY);
    newRed = minR + ((maxR - minR)/(maxY - minY))*(splitY - minY);
    newGreen = minG + ((maxG - minG)/(maxY - minY))*(splitY - minY);
    newBlue = minB + ((maxB - minB)/(maxY - minY))*(splitY - minY);
    newNormX = minNormX + ((maxNormX - minNormX)/(maxY - minY))*(splitY - minY);
    newNormY = minNormY + ((maxNormY - minNormY)/(maxY - minY))*(splitY - minY);
    newNormZ = minNormZ + ((maxNormZ - minNormZ)/(maxY - minY))*(splitY - minY);
    Triangle upwardTriangle,downwardTriangle;
    upwardTriangle.normals[0][0] = maxNormX;
    upwardTriangle.normals[0][1] = maxNormY;
    upwardTriangle.normals[0][2] = maxNormZ;
    downwardTriangle.normals[0][0] = minNormX;
    downwardTriangle.normals[0][1] = minNormY;
    downwardTriangle.normals[0][2] = minNormZ;

    upwardTriangle.normals[1][0] = downwardTriangle.normals[1][0] = splitNormX;
    upwardTriangle.normals[1][1] = downwardTriangle.normals[1][1] = splitNormY;
    upwardTriangle.normals[1][2] = downwardTriangle.normals[1][2] = splitNormZ;

    upwardTriangle.normals[2][0] = downwardTriangle.normals[2][0] = newNormX;
    upwardTriangle.normals[2][1] = downwardTriangle.normals[2][1] = newNormY;
    upwardTriangle.normals[2][2] = downwardTriangle.normals[2][2] = newNormZ;

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

int main(){
    int width = 1000, height = 1000;
    vtkImageData *image = NewImage(width, height);
    unsigned char *buffer =
      (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = width*height;
    double zvalues[npixels];

    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.buffer = buffer;
    screen.width = width;
    screen.height = height;

    Camera camera;
    for (int i = 0; i < 1000; i+=1000){
        //initialize values
        std::fill_n(zvalues, npixels, -1.0);
        std::fill_n(buffer, npixels*3, 0);

        camera = GetCamera(i, 1000);

        camera.cameraTransform();
        camera.viewTransform();
        camera.deviceTransform(screen);
        Matrix matrix = getObjectMatrix(camera,screen);

        //transform triangles based off of camera
        for (int j = 0; j < triangles.size(); j++){
            Triangle triangle = TriangleTransform(triangles[j], matrix);
            for (int n = 0; n <= 2; n++){ angleView[n] = camera.v3[n]; }
            placeTrianglePixels(triangle,buffer,zvalues,width,height);
        }
        //name files and write the image
        char filename[32];
        sprintf(filename, "frame%03d", i);
        WriteImage(image, filename);
    }
}
