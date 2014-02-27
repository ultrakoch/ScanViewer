/**
* @file   ${fileShow.cpp}
* @author Rainer Koch
* @date   ${2014/18/02}
*
* Load different file clouds to filter reflections
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <memory>
#include <stdlib.h>
#include <math.h>
//#include "VisualObjects.h"

#include "fileShow.h"
#include "obgraphic/Obvious3D.h"
#include <vector>

#define _USE_MATH_DEFINES

using namespace std;
using namespace obvious;

/* Default parameters for file*/
   char* file;
   string filetype;
   unsigned int file_nr = 0;
   bool _save = true;

 /* Default parameters for viewer */
   std::vector<long> dist;
   double* intens;
   VtkCloud* cloud_scan;
   VtkCloud* cloud_show;
   VtkCloud* cloud_mirror;
   VtkCloud* cloud_sensor;

   int reducePoints=5;

   int cloudsize_scan;
   int cloudsize_show;

   int raysPscan;
   int sensorsize = 5;  // only odd numbers
   double* data_scan;
   double* data_show;
   double* data_mirror;
   double* data_sensor;
   double* intensity_scan;

   double* sensorPos;
   unsigned char* colors_scan;
   unsigned char* colors_show;
   unsigned char* colors_mirror;
   unsigned char* colors_sensor;
   unsigned char intensColor;

   Obvious3D* viewer3D = new Obvious3D();
   bool _pause = false;
   bool _filterSwitch = false;
   bool _pointReducer = false;
   bool _pointRIntensityShow = false;

   int id = 0;            // hokuyo measurement id

   // Colorlevels
   int cLevel1 = 3000;
   int cLevel2 = 3500;

   // boundry box for measurements
   // values are limits for pos. and neg. distance
   double max_x = 200.0;
   double max_y = 200.0;
   double max_z = 200.0;

   // Filter variables
   double lineColor[3];
   double mirrorPlaneColor[3];
   unsigned int mirrorPlaneSize_X = 1;
   unsigned int mirrorPlaneSize_Y = 1;
   double* minIntensCoord;
   unsigned char* minIntensColor;
   double* maxIntensCoord;
   unsigned char* maxIntensColor;
   double min_intensity;
   double max_intensity;
   double** maxIntensityPos;


   // Test variables
   unsigned int testFktCount = 0;

   bool test = 0;
   bool test1 = 0;
   int max_tmp;
   int min_tmp;
   unsigned int tmp = 0;



void testfkt()
{
  cout << " was here " << testFktCount << endl;
  testFktCount++;
}

void init(int argc, char* argv[])
{
  if(argc<2)
  {
     cout << "usage: " << argv[0] << "[filename] [txt_xy, txt_xyi, txt_3xy, txt_3xyi, rxp ]" << endl;
  }
  else
  {
   file = argv[1];
   filetype = argv[2];
  }

  // Variables:
  lineColor[0] = 100;
  lineColor[1] = 100;
  lineColor[2] = 100;

  mirrorPlaneColor[0] = 50;
  mirrorPlaneColor[1] = 50;
  mirrorPlaneColor[2] = 0;

  minIntensCoord   = new double[3];
  minIntensColor   = new unsigned char[3];
  maxIntensCoord   = new double[3];
  maxIntensColor   = new unsigned char[3];
  double** sensorPos = new double*[1];
  sensorPos[0]= new double[3];
  sensorPos[0][0] = 0.0;
  sensorPos[0][1] = 0.0;
  sensorPos[0][2] = 0.0;

  maxIntensityPos = new double* [1];
  maxIntensityPos[0]= new double[3];
}

void initShowCloud()
{
  if(_pointReducer)
  {
    cloudsize_show = cloudsize_scan - (cloudsize_scan /  reducePoints);

    cout << "cloud_scan: " << cloudsize_scan << endl;
    cout << "cloud_show: " << cloudsize_show << endl;
    cout << "reduce every: " << reducePoints << endl;
    cout << "Point reducer: " << _pointReducer << endl;

  data_show   = new double[cloudsize_show * 3];
  colors_show = new unsigned char[cloudsize_show * 3];
  }
  else
  {
    cloudsize_show = cloudsize_scan;
    data_show   = new double[cloudsize_show * 3];
    colors_show = new unsigned char[cloudsize_show * 3];
  }
}
void createShowCloud()
{
  int i=0;
  int f=0;
  int n=0;
  // Point reducer to show less points than necessary
  if(_pointReducer)
  {
    while(i < cloudsize_scan)
    {
      data_show[3*f] = data_scan[3*i];
      data_show[3*f+1] = data_scan[3*i+1];
      data_show[3*f+2] = data_scan[3*i+2];

      colors_show[3*f] = colors_scan[3*i];
      colors_show[3*f+1] = colors_scan[3*i+2];
      colors_show[3*f+2] = colors_scan[3*i+2];

      f++;
      i = i+1;
      if(n==reducePoints)
        n = 0;
      else
      {
        n++;
        i = i+1;
      }
    }
  }
  else
  {
    data_show = data_scan;
    colors_show = colors_scan;
  }

  cloud_show->setCoords(data_show, cloudsize_show, 3);
  cloud_show->setColors(colors_show, cloudsize_show, 3);
}

void createCube(int cubesize)
{
  //VtkCloud* cloud_cube;
  //cloud_cube = new VtkCloud();
  double* data_cube;
  data_cube   = new double[3*cubesize*cubesize*cubesize];
  unsigned char* color_cube;
  color_cube = new unsigned char[3*cubesize*cubesize*cubesize];


  int n = 0;
  for (int i=0; i<=cubesize-1; i++)
  {
    for (int f=0; f<=cubesize-1; f++)
    {
      for (int l=0; l<=cubesize-1; l++)
      {
        data_cube[3*n] = (cubesize/2)-i;
        data_cube[3*n+1] = (cubesize/2)-f;
        data_cube[3*n+2] = (cubesize/2)-l;

        color_cube[3*n] = 255;
        color_cube[3*n+1] = 215;
        color_cube[3*n+2] = 0;
     //   cout << data_cube[3*n] << " " <<  data_cube[3*n+1] << " " << data_cube[3*n+2] << endl;
        n++;
      }
    }

  }

  cloud_sensor->setCoords(data_cube, (cubesize*cubesize*cubesize), 3);
  cloud_sensor->setColors(color_cube, (cubesize*cubesize*cubesize), 3);
}

void deltectMirroredPoints(double* cloud, int size, double origin[3], double axis1[3], double axis2[3])
{
  // line[] = origin[] + lamda * (origin[]-axis1[])
  // => lamda = (line[] - origin[] ) / (origin[]-axis1[])
  double lamda;
  double n_mirror[3];
  double v_axis1[3];
  double v_axis2[3];
  double line[3];

  // calculate normal vector of mirror
  v_axis1[0] = axis1[0] - origin[0];
  v_axis1[1] = axis1[1] - origin[1];
  v_axis1[2] = axis1[2] - origin[2];

  v_axis2[0] = axis2[0] - origin[0];
  v_axis2[1] = axis2[1] - origin[1];
  v_axis2[2] = axis2[2] - origin[2];

  n_mirror[0] = v_axis1[1] * v_axis2[2] - v_axis1[2] * v_axis2[1];
  n_mirror[1] = v_axis1[2] * v_axis2[0] - v_axis1[0] * v_axis2[2];
  n_mirror[2] = v_axis1[0] * v_axis2[1] - v_axis1[1] * v_axis2[0];


  for(int i=0; i < size; i++)
  {
    // calculate lamda / point on the mirror line at x value of cloud
    if(n_mirror[1] != 0) // otherwise mirror is parallel to x axis
    {
      lamda = (cloud[3*i] - origin[0]) / (origin[0] - axis1[0]);
      // calculate y value on the line
      line[1] = origin[1] + (lamda * (origin[1] - axis1[1]));

// TODO: Expand to check also z-value => check for plane
// TODO: Expand to check only for lamda < 1 => origin, axis1 and axis2 are exactly corner points of the mirror! => mirror is not a endless plane anymore, but a fixed size plane
      // Check points on location
      if(cloud[3*i+1] == line[1]) // Point on mirror line
      {
        data_mirror[3*i] = cloud[3*i];
        data_mirror[3*i+1] = cloud[3*i+1];
        data_mirror[3*i+2] = cloud[3*i+2];

        colors_mirror[3*i]     = 0;                         // r
        colors_mirror[3*i +1]   = 0;                     // g
        colors_mirror[3*i +2]   = 250;                         // b
      }
      else if((cloud[3*i+1] > line[1]) && (n_mirror[1] > 0) or
              (cloud[3*i+1] < line[1]) && (n_mirror[1] < 0))      // Point behind mirror plane
      {
        data_mirror[3*i] = cloud[3*i];
        data_mirror[3*i+1] = cloud[3*i+1];
        data_mirror[3*i+2] = cloud[3*i+2];

        colors_mirror[3*i]     = 250;                         // r
        colors_mirror[3*i +1]   = 0;                     // g
        colors_mirror[3*i +2]   = 0;                         // b
      }
      else  // regular point
      {
        data_mirror[3*i] = 0;
        data_mirror[3*i+1] = 0;
        data_mirror[3*i+2] = 0;

        colors_mirror[3*i]     = 0;                         // r
        colors_mirror[3*i +1]   = 0;                     // g
        colors_mirror[3*i +2]   = 0;                         // b
       }
    }
    else if(n_mirror[1] == 0)
    {
      line[0] = origin[0];
      // Check points on location
      if(cloud[3*i] == line[0]) // Point on mirror line
      {
        data_mirror[3*i] = cloud[3*i];
        data_mirror[3*i+1] = cloud[3*i+1];
        data_mirror[3*i+2] = cloud[3*i+2];

        colors_mirror[3*i]     = 0;                         // r
        colors_mirror[3*i +1]   = 0;                     // g
        colors_mirror[3*i +2]   = 250;                         // b
      }
      else if((cloud[3*i] > line[0]) && (n_mirror[0] > 0) or
          (cloud[3*i] < line[0]) && (n_mirror[0] < 0))
      {
        data_mirror[3*i] = cloud[3*i];
        data_mirror[3*i+1] = cloud[3*i+1];
        data_mirror[3*i+2] = cloud[3*i+2];

        colors_mirror[3*i]     = 250;                         // r
        colors_mirror[3*i +1]   = 0;                     // g
        colors_mirror[3*i +2]   = 0;                         // b
      }
      else  // regular point
      {
        data_mirror[3*i] = 0;
        data_mirror[3*i+1] = 0;
        data_mirror[3*i+2] = 0;

        colors_mirror[3*i]     = 0;                         // r
        colors_mirror[3*i +1]   = 0;                     // g
        colors_mirror[3*i +2]   = 0;                         // b
       }
    }
    else
    {
      cout << "Mirro plane not possible 1" << endl;
    }
   }
}


void filter(double* distance, unsigned char* colors, double* intensity, int cloudsize)
{
  float distance_refl = 0.0;
  float angle_refl = 0.0;
  float angle_refl_g = 0.0;

  float distance_orig = 0.0;
  float angle_orig = 0.0;
  float angle_orig_g = 0.0;

  unsigned int id_orig = 0;

// Find maximum and minimum intensity
  min_intensity = intensity[0];
  minIntensCoord[0] = distance[0];
  minIntensCoord[0+1] = distance[1];
  minIntensCoord[0+2] = distance[2];
  minIntensColor[0] = colors[0];
  minIntensColor[0+1] = colors[1];
  minIntensColor[0+2] = colors[2];
  max_intensity = intensity[0];
  maxIntensCoord[0] = distance[0];
  maxIntensCoord[0+1] = distance[1];
  maxIntensCoord[0+2] = distance[2];
  maxIntensColor[0] = colors[0];
  maxIntensColor[0+1] = colors[1];
  maxIntensColor[0+2] = colors[2];

  for(int i=1; i<cloudsize; i++)
  {
    if(intensity[i] < min_intensity)
    {
      min_intensity = intensity[i];
      minIntensCoord[0] = distance[3*i];
      minIntensCoord[0+1] = distance[3*i+1];
      minIntensCoord[0+2] = distance[3*i+2];
      minIntensColor[0] = colors[3*i];
      minIntensColor[0+1] = colors[3*i+1];
      minIntensColor[0+2] = colors[3*i+2];
    }
    if(intensity[i] > max_intensity)
    {
      max_intensity = intensity[i];
      maxIntensCoord[0] = distance[3*i];
      maxIntensCoord[0+1] = distance[3*i+1];
      maxIntensCoord[0+2] = distance[3*i+2];
      maxIntensColor[0] = colors[3*i];
      maxIntensColor[0+1] = colors[3*i+1];
      maxIntensColor[0+2] = colors[3*i+2];
    }
  }

  //cout << "Min: " << min_intensity << " Koordinaten: " << minIntensCoord[0] << " / " << minIntensCoord[1] << " / " << minIntensCoord[2] << endl;
  //cout << "Max: " << max_intensity << " Koordinaten: " << maxIntensCoord[0] << " / " << maxIntensCoord[1] << " / " << maxIntensCoord[2] << endl;


// add line to total reflection point
  maxIntensityPos[0][0] = maxIntensCoord[0];
  maxIntensityPos[0][1] = maxIntensCoord[1];
  maxIntensityPos[0][2] = maxIntensCoord[2];

  viewer3D->addLines(&sensorPos, maxIntensityPos, 1, lineColor);

// Mirror Plane
  /* Points to build up mirror plane from max Intensity Point */
//  double point1[3];
//  point1[0] = maxIntensCoord[0]-maxIntensCoord[1];
//  point1[1] = maxIntensCoord[0]+maxIntensCoord[1];
//  point1[2] = 0;
//  double point2[3];
//  point2[0] = maxIntensCoord[0];
//  point2[1] = maxIntensCoord[1];
//  point2[2] = 50;
//  // move origin of mirror
//  double origin[3];
//  origin[0] = maxIntensCoord[0];
//  origin[1] = maxIntensCoord[1];
//  origin[2] = maxIntensCoord[2];

  /* Points to build up mirror plane with max Intensity Point at center */
  double origin[3];
  origin[0] = maxIntensCoord[0]+maxIntensCoord[1];
  origin[1] = -maxIntensCoord[0]-maxIntensCoord[1];
  origin[2] = -50;
  double point1[3];
  point1[0] = maxIntensCoord[0]-maxIntensCoord[1];
  point1[1] = maxIntensCoord[0]+maxIntensCoord[1];
  point1[2] = -50;
  double point2[3];
  point2[0] = origin[0];
  point2[1] = origin[1];
  point2[2] = 50;

//  cout << "Origin/Max: " << " Koordinaten: " << maxIntensCoord[0] << " / " << maxIntensCoord[1] << " / " << maxIntensCoord[2] << endl;
//  cout << "Point 1: " << " Koordinaten: " << point1[0] << " / " << point1[1] << " / " << point1[2] << endl;
//  cout << "Point 2: " << " Koordinaten: " << point2[0] << " / " << point2[1] << " / " << point2[2] << endl;

  viewer3D->addPlane(origin, point1, point2, mirrorPlaneSize_X, mirrorPlaneSize_Y, mirrorPlaneColor);

// find points behind mirror plane
   deltectMirroredPoints(data_scan, cloudsize_scan, origin, point1, point2);

  /* mark points behind mirror plan */
  cloud_mirror->setCoords(data_mirror, cloudsize, 3);
  cloud_mirror->setColors(colors_mirror, cloudsize, 3);

  }
void boundryBox()
{

  for(int i=0; i<cloudsize_scan; i++)
  {
    //cout << data_scan[3*i] << " " << data_scan[3*i+1] << " " << data_scan[3*i+2] << endl;
    if((abs(data_scan[3*i]) >= max_x) or (abs(data_scan[3*i+1]) >= max_y) or (abs(data_scan[3*i+2]) >= max_z))
    {
      //cout << "was here " << endl;
      data_scan[3*i]       = 0.0;
      data_scan[3*i + 1]   = 0.0;
      data_scan[3*i + 2]   = 0.0;

      colors_scan[3*i]     = 0;                     // r
      colors_scan[3*i+1]   = 0;                     // g
      colors_scan[3*i+2]   = 0;                     // b

      intensity_scan[i] = 0;
    }
    //cout << data_scan[3*i] << " " << data_scan[3*i+1] << " " << data_scan[3*i+2] << endl;
    i++;
  }

}

int load_xy_file(char* filename)
{
  // Read txt file with x, y
      float x, y, z, degree, intensity;
      string leseString;
      ifstream file(filename);
      int echoCount = 0;
      unsigned int linesFile = 0;
      unsigned int id = 0;

      /* read single distance */
      if(!filename)
        cerr << "Eingabe-Datei kann nicht geöffnet werden\n";
      else
      {
        for(leseString; getline(file, leseString);)
        {
          linesFile++;
        }
        id = linesFile -1;
        cout << "lines: " << linesFile << endl;
        linesFile = 0;
        data_scan   = new double[id * 3]; // linesFile - header
        colors_scan = new unsigned char[id * 3];
        std::ifstream file(filename);
        // print file header and remove it
        getline(file, leseString);
        cout << "Fileheader:" << leseString << endl;
        for(leseString; getline(file, leseString);)
        {
          istringstream(leseString,ios_base::in) >> x >> y;
          // Build vtkCloud
          data_scan[3*linesFile]       = x;  // x
          data_scan[3*linesFile + 1]   = y;  // y
          data_scan[3*linesFile + 2]   = 0;                                  // z

          colors_scan[3*linesFile]     = 200;                                // r
          colors_scan[3*linesFile+1]   = 200;                                // g
          colors_scan[3*linesFile+2]   = 200;                                // b
         // cout << "data id: " << linesFile << " x: " << data[3*linesFile] << "y: " << data[3*linesFile+1] << endl;
          linesFile++;
         }
      }
      cout << "File " << filename << " with x,y loaded" << endl;
        return id;
}
int load_xyi_file(char* filename)
{
  // Read txt file with x, y
      float x, y, z, degree, intensity;
      string leseString;
      ifstream file(filename);
      int echoCount = 0;
      unsigned int linesFile = 0;
      unsigned int id = 0;

      /* read single distance and intensity */
      if(!filename)
        cerr << "Eingabe-Datei kann nicht geöffnet werden\n";
      else
      {
        for(leseString; getline(file, leseString);)
        {
          linesFile++;
        }
        id = linesFile -1;
        cout << "lines: " << linesFile << endl;
        linesFile = 0;
        data_scan   = new double[id * 3]; // linesFile - header
        intensity_scan   = new double[id]; // linesFile - header
        colors_scan = new unsigned char[id * 3];
        std::ifstream file(filename);
        // print file header and remove it
        getline(file, leseString);
        cout << "Fileheader:" << leseString << endl;
        for(leseString; getline(file, leseString);)
        {
          istringstream(leseString,ios_base::in) >> x >> y >> intensity;

          // Build vtkCloud
          data_scan[3*linesFile]       = x;  // x
          data_scan[3*linesFile + 1]   = y;  // y
          data_scan[3*linesFile + 2]   = 0;                                  // z

          intensity_scan[linesFile] = intensity;

          colors_scan[3*linesFile]     = 100;                         // r
          colors_scan[3*linesFile+1]   = 100;                     // g
          colors_scan[3*linesFile+2]   = 100;                         // b

        //  cout << "data_scan id: " << linesFile << " x: " << data_scan[3*linesFile] << " y: " << data_scan[3*linesFile+1] << " int: " << intensity_scan[linesFile] << endl;
          linesFile++;
         }
      }
      cout << "File " << filename << " with x,y,intensity loaded" << endl;
     return id;
}
int load_3xy_file(char* filename)
{
  // Read txt file with x, y
      float x, y, z, degree, intensity;
      string leseString;
      ifstream file(filename);
      int echoCount = 0;
      unsigned int linesFile = 0;
      unsigned int id = 0;

      /* read multi distance */
      if(!filename)
        cerr << "Eingabe-Datei kann nicht geöffnet werden\n";
      else
      {
        for(leseString; getline(file, leseString);)
        {
          linesFile++;
        }
        id = linesFile -1;
        linesFile = 0;
        //cout << "lines: " << linesFile << endl;
        data_scan   = new double[id * 3]; // linesFile - header
        colors_scan = new unsigned char[id * 3];
        std::ifstream file(filename);
        // print file header and remove it
        getline(file, leseString);
        cout << "Fileheader:" << leseString << endl;
        for(leseString; getline(file, leseString);)
        {
          istringstream(leseString,ios_base::in) >> x >> y;

          // Build vtkCloud
          data_scan[3*linesFile]       = x;  // x
          data_scan[3*linesFile + 1]   = y;  // y

          if(echoCount==0)
          {
            // First echo
            data_scan[3*linesFile + 2]   = 0;                                  // z

            colors_scan[3*linesFile]     = 200;                                // r
            colors_scan[3*linesFile+1]   = 0;                                // g
            colors_scan[3*linesFile+2]   = 0;                                // b
          }
          else if(echoCount==1)
          {
            // Second echo
            data_scan[3*linesFile + 2]   = 25;                                  // z

            colors_scan[3*linesFile]     = 0;                                // r
            colors_scan[3*linesFile+1]   = 200;                                // g
            colors_scan[3*linesFile+2]   = 0;                                // b
          }
          else if(echoCount==2)
          {
            // Third echo
            data_scan[3*linesFile + 2]   = 50;                                  // z

            colors_scan[3*linesFile]     = 0;                                // r
            colors_scan[3*linesFile+1]   = 0;                                // g
            colors_scan[3*linesFile+2]   = 200;                                // b
          }
          if(echoCount < 2)
            echoCount++;
          else
            echoCount = 0;
          //cout << "data_scan id: " << linesFile << " x: " << data_scan[3*linesFile] << "y: " << data_scan[3*linesFile+1] << endl;
          linesFile++;
         }
      }

        return id;
}
int load_3xyi_file(char* filename)
{
  // Read txt file with x, y
      float x, y, z, degree, intensity;
      string leseString;
      ifstream file(filename);
      int echoCount = 0;
      max_tmp = 32767;
      unsigned int linesFile = 0;
      unsigned int id = 0;

      /* read multi distance and intensity */
      if(!filename)
        cerr << "Eingabe-Datei kann nicht geöffnet werden\n";
      else
      {
        for(leseString; getline(file, leseString);)
        {
          linesFile++;
        }
        id = linesFile -1;
        linesFile = 0;
        //cout << "lines: " << linesFile << endl;
        data_scan   = new double[id * 3]; // linesFile - header
        colors_scan = new unsigned char[id * 3];
        std::ifstream file(filename);
        // print file header and remove it
        getline(file, leseString);
        cout << "Fileheader:" << leseString << endl;
        for(leseString; getline(file, leseString);)
        {
          istringstream(leseString,ios_base::in) >> x >> y >> intensity;

          // Build vtkCloud
          data_scan[3*linesFile]       = x;  // x
          data_scan[3*linesFile + 1]   = y;  // y
          intensity_scan[linesFile] = intensity;

          intensColor = intensity/max_tmp*256;

          if(echoCount==0)
          {
            // First echo
            data_scan[3*linesFile + 2]   = 0;                                  // z

            colors_scan[3*linesFile]     = intensColor;                         // r
            colors_scan[3*linesFile+1]   = 100;                     // g
            colors_scan[3*linesFile+2]   = 256-intensColor;                         // b
          }
          else if(echoCount==1)
          {
            // Second echo
            data_scan[3*linesFile + 2]   = 25;                                  // z

            colors_scan[3*linesFile+3]   = 256-intensColor;                         // r
            colors_scan[3*linesFile+4]   = intensColor;                     // g
            colors_scan[3*linesFile+5]   = 100;                         // b
          }
          else if(echoCount==2)
          {
            // Third echo
            data_scan[3*linesFile + 2]   = 50;                                  // z

            colors_scan[3*linesFile+6]   = 100;                         // r
            colors_scan[3*linesFile+7]   = 256-intensColor;                     // g
            colors_scan[3*linesFile+8]   = intensColor;                         // b
          }
          if(echoCount < 2)
            echoCount++;
          else
            echoCount = 0;
          //cout << "data_scan id: " << linesFile << " x: " << data_scan[3*linesFile] << "y: " << data_scan[3*linesFile+1] << endl;
          linesFile++;
         }
      }

        return id;
}
int load_rxp_file(char* filename)
{
  // Read rxp file
      float x, y, z, degree, intensity, reflectance, deviation;
      string leseString;
      ifstream file(filename);
      int echoCount = 0;
      max_tmp = 32767;
      unsigned int linesFile = 0;
      unsigned int id = 0;

      if(!filename)
         cerr << "Eingabe-Datei kann nicht geöffnet werden\n";
       else
       {
         for(leseString; getline(file, leseString);)
         {
           linesFile++;
         }

         id = linesFile -2; // Header have two lines
         //cout << "lines: " << linesFile << endl;
         linesFile = 0;
         data_scan   = new double[id * 3]; // linesFile - header
         intensity_scan   = new double[id]; // linesFile - header
         colors_scan = new unsigned char[id * 3];
         std::ifstream file(filename);

         // print file header and remove it
         getline(file, leseString);
         cout << "Fileheader:" << leseString << endl;
         getline(file, leseString);
         cout << "Fileheader:" << leseString << endl;
         for(leseString; getline(file, leseString);)
         {
           istringstream(leseString,ios_base::in) >> x >> y >> z >> reflectance >> intensity >> deviation;

           // Build vtkCloud
           data_scan[3*linesFile]       = x;
           data_scan[3*linesFile + 1]   = y;
           data_scan[3*linesFile + 2]   = z;

           intensity_scan[linesFile] = reflectance;

// TODO: reflectance and deviation not used yet

           if(intensity<= cLevel1)
           {
             colors_scan[3*linesFile]     = 100;                         // r
             colors_scan[3*linesFile+1]   = 0;                     // g
             colors_scan[3*linesFile+2]   = 0;                         // b
           }
           else if ((intensity > cLevel1) && (intensity<= cLevel2))
           {
             colors_scan[3*linesFile]     = 0;                         // r
             colors_scan[3*linesFile+1]   = 100;                     // g
             colors_scan[3*linesFile+2]   = 0;                         // b
           }
           else if (intensity > cLevel1)
           {
             colors_scan[3*linesFile]     = 100;                         // r
             colors_scan[3*linesFile+1]   = 100;                     // g
             colors_scan[3*linesFile+2]   = 100;                         // b

             int abs_id = abs(id-512);
             colors_scan[3*abs_id ]     = 100;                         // r
             colors_scan[3*abs_id +1]   = 100;                     // g
             colors_scan[3*abs_id +2]   = 100;                         // b

           }
           linesFile++;
          }
       }
      cout << "File " << filename << " as rxp with x,y,intensity loaded" << endl;
      return id;
}

class vtkTimerCallback : public vtkCommand
{
public:
  static vtkTimerCallback *New()
  {
    vtkTimerCallback *cb = new vtkTimerCallback;
    return cb;
  }

  virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long eventId,  void *vtkNotUsed(callData))
  {

/* Filter */
      if(_filterSwitch == 1)
        filter(data_scan, colors_scan, intensity_scan, cloudsize_scan);

/* Set new cloud*/
      createShowCloud();
      viewer3D->update();

    }

private:

};
void cbPointIncrease()
{
  if(reducePoints>5)
  {
    reducePoints = reducePoints-1;
    cout << "Increasing " << reducePoints << " Points!" << endl;
    initShowCloud();
  }
  else
  {
    cout << "No more increasing possible (" << reducePoints << ")" << endl;
  }
}
void cbPointReduce()
{
  if(reducePoints<100)
  {
    reducePoints = reducePoints+1;
    cout << "Reduceing " << reducePoints << " Points!" << endl;
    initShowCloud();
  }
  else
  {
    cout << "No more reduceing possible (" << reducePoints << ")" << endl;

  }
}
void cbPointReducerEnable()
{
  _pointReducer = !_pointReducer;
  if(_pointReducer)
    cout << "Point reduction on!" << endl;
  if(!_pointReducer)
    cout << "Point reduction off!" << endl;
  initShowCloud();
}
void cbIntensityColor()
{
  _pointRIntensityShow = !_pointRIntensityShow;
  if(_pointRIntensityShow)
  {
    for(int i=0; i<cloudsize_scan; i++)
    {
      if(intensity_scan[i]<= cLevel1)
      {
        colors_scan[3*i]     = 100;                         // r
        colors_scan[3*i+1]   = 100;                     // g
        colors_scan[3*i+2]   = 100;                         // b
      }
      else if ((intensity_scan[i] > cLevel1) && (intensity_scan[i]<= cLevel2))
      {
        colors_scan[3*i]     = 0;                         // r
        colors_scan[3*i+1]   = 200;                     // g
        colors_scan[3*i+2]   = 0;                         // b
      }
      else if (intensity_scan[i] > cLevel2)
      {
        colors_scan[3*i]     = 200;                         // r
        colors_scan[3*i+1]   = 0;                     // g
        colors_scan[3*i+2]   = 0;                         // b
       }
    }
  }
  else
  {
    for(int i=0; i<cloudsize_scan; i++)
      {
        colors_scan[3*i]     = 100;                         // r
        colors_scan[3*i+1]   = 100;                     // g
        colors_scan[3*i+2]   = 100;                         // b
      }
  }
}

int main(int argc, char* argv[])
{

/* Init */
  // Load parameters
  init(argc, argv);
  cout << "init finished" << endl;

  // Initialization of 3D viewer
  cloud_show = new VtkCloud();
  cloud_mirror = new VtkCloud();
  cloud_sensor = new VtkCloud();
  // define update timer
   vtkSmartPointer<vtkTimerCallback> cb =  vtkSmartPointer<vtkTimerCallback>::New();
   vtkSmartPointer<vtkRenderWindowInteractor> interactor = viewer3D->getWindowInteractor();
   interactor->AddObserver(vtkCommand::TimerEvent, cb);
   interactor->CreateRepeatingTimer(30);
   // define keyboard callbacks
   viewer3D->registerFlipVariable("f", &_filterSwitch);
   viewer3D->registerKeyboardCallback("1", cbPointReducerEnable); // enable and disable point reduction
   viewer3D->registerKeyboardCallback("d", cbPointReduce); // reduce amount of shown points
   viewer3D->registerKeyboardCallback("i", cbPointIncrease); // increase amount of shown points
   viewer3D->registerKeyboardCallback("c", cbIntensityColor); // color the data points depending on their intensity

   cout << "VTK vierwer init" << endl;


/* Load data */
   if(filetype.compare("rxp")==0)
      cloudsize_scan = load_rxp_file(file);
   else if(filetype.compare("txt_xy")==0)
     cloudsize_scan = load_xy_file(file);
   else if(filetype.compare("txt_xyi")==0)
     cloudsize_scan = load_xyi_file(file);
   else if(filetype.compare("txt_3xy")==0)
     cloudsize_scan = load_3xy_file(file);
   else if(filetype.compare("txt_3xyi")==0)
     cloudsize_scan = load_3xyi_file(file);
   else
     cout << "Can not load filetype " << filetype << endl;

/* Prefilter data */
//TODO: Not working yet
   //boundryBox();


/* Create objects */
   initShowCloud();
   createShowCloud();
   // Create mirror plane
   // size = data_size_Scanner * 3_koordinates
   data_mirror   = new double[cloudsize_scan * 3];
   colors_mirror = new unsigned char[cloudsize_scan * 3];
   // Show sensor cube
   sensorPos =new double[2];
   sensorPos[0] = 0.0;
   sensorPos[1] = 0.0;
   sensorPos[2] = 0.0;
   createCube(sensorsize);
   viewer3D->addCloud(cloud_sensor);
   cout << "Cube created" << endl;
   // Show axes
   viewer3D->showAxes(true);
   cout << "Axes created" << endl;

/* Load and start vierwer */
   viewer3D->addCloud(cloud_show);
   viewer3D->addCloud(cloud_mirror);
   viewer3D->startRendering();


/* Shutdown */
   delete viewer3D;

   delete [] data_scan;
   delete [] data_show;
   delete [] data_mirror;
   delete [] data_sensor;

   delete [] colors_show;
   delete [] colors_scan;
   delete [] colors_mirror;
   delete [] colors_sensor;

   delete cloud_show;
   delete cloud_mirror;
   delete cloud_sensor;
   return 0;

}
