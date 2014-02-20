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

    // 0 = Laser, 1 = File
   int sourcetype = 0;

/* Default parameters for file*/
   char* file;
   string filetype;
   unsigned int file_nr = 0;
   bool _save = true;

 /* Default parameters for viewer */
   std::vector<long> dist;
   float* intens;
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
   bool _pointReducer = true;

   int id = 0;            // hokuyo measurement id

   // Testvariables
   bool test = 0;
   bool test1 = 0;
   unsigned int max_tmp = 0;
   unsigned int min_tmp = 65535;
   unsigned int tmp = 0;

   // Colorlevels
   int cLevel1 = 200;
   int cLevel2 = 400;

   // boundry box for measurements
   // values are limits for pos. and neg. distance
   double max_x = 100.0;
   double max_y = 100.0;
   double max_z = 100.0;


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
}

void initShowCloud()
{
  if(_pointReducer)
    cloudsize_show = cloudsize_scan - (cloudsize_scan /  reducePoints);
  else
    cloudsize_show = cloudsize_scan;

  cout << "cloud_scan: " << cloudsize_scan << endl;
  cout << "cloud_show: " << cloudsize_show << endl;
  cout << "reduce every: " << reducePoints << endl;
  cout << "Point reducer: " << _pointReducer << endl;

  data_show   = new double[cloudsize_show * 3];
  colors_show = new unsigned char[cloudsize_show * 3];
}
void createShowCloud()
{
  int i=0;
  int f=0;
  int n=0;
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

void filter(double* distance, unsigned char* colors, float* intensity)
{

  float distance_refl = 0.0;
  float angle_refl = 0.0;
  float angle_refl_g = 0.0;

  float distance_orig = 0.0;
  float angle_orig = 0.0;
  float angle_orig_g = 0.0;

  unsigned int id_orig = 0;


    for(int i=0; i<=cloudsize_scan; i++)
    {
      // threshold for intensity
      if(intensity[i] > 2700)
      {
        // recalculate original point id
        id_orig = abs(cloudsize_scan-i);
//           cout << "id_refl: " << i << " id_orig: " << id_orig << endl;

        // calculate vector
        //TODO: later expand with z => sqrt(x2+y2+z2)
        distance_refl = sqrt(((distance[3*i])*(distance[3*i])) + (distance[3*i+1]*distance[3*i+1]));
        angle_refl_g = asin(distance[3*i]/distance_refl*1.0)/M_PI*360;
        angle_refl = asin(distance[3*i]/distance_refl*1.0);

        distance_orig = sqrt(((distance[3*id_orig])*(distance[3*id_orig])) + (distance[3*id_orig+1]*distance[3*id_orig+1]));
        angle_orig_g = abs(360.0-angle_refl);
        angle_refl = M_PI-asin(distance[3*i]/distance_refl*1.0);


        // mark reflected point
        colors[3*i]     = 100;                         // r
        colors[3*i +1]   = 0;                     // g
        colors[3*i +2]   = 0;                         // b

//        // mark original point
//        colors_scan[3*id_orig]     = 100;                         // r
//        colors_scan[3*id_orig +1]   = 0;                     // g
//        colors_scan[3*id_orig +2]   = 100;                         // b

        // calculate mirror plane
        data_mirror[3*i] = abs(distance_orig-distance_refl)*sin(angle_refl);
        data_mirror[3*i+1] = abs(distance_orig-distance_refl)*cos(angle_refl);
        data_mirror[3*i+2] = 0;

        colors_mirror[3*i]     = 50;                         // r
        colors_mirror[3*i +1]   = 200;                     // g
        colors_mirror[3*i +2]   = 50;                         // b

        if(test==0)
        {
          cout << " dist ref: " << distance_refl << " dist orig: " << distance_orig << endl;
          cout << " id_refl: " << i << " id_orig: " << id_orig << endl;
          cout << " angl_refl: " << angle_refl << " angl_orig: " << angle_orig << endl;
          cout << " angl_refl: " << angle_refl_g << "° angl_orig: " << angle_orig_g << "°" << endl;
          cout << " x_refl: " << distance[3*i] << " y_ref: " << distance[3*i+1] << endl;
          //cout << " x_orig: " << distance[3*id_orig] << " y_orig: " << distance[3*id_orig+1] << endl;
          cout << " x_mirr: " << data_mirror[3*i] << " y_mirr: " << data_mirror[3*i+1] << endl << endl << endl;
        }
      }

    test = 1;

    /* Set recalculated mirro plane */
    cloud_mirror->setCoords(data_mirror, cloudsize_scan, 3);
    cloud_mirror->setColors(colors_mirror, cloudsize_scan, 3);

  }

  }
void boundryBox()
{
  for(int i=0; i<cloudsize_scan; i++)
  {
    cout << data_scan[3*i] << " " << data_scan[3*i+1] << " " << data_scan[3*i+2] << endl;
    if((abs(data_scan[3*i]) >= max_x) or (abs(data_scan[3*i+1]) >= max_y) or (abs(data_scan[3*i+2]) >= max_z))
    {
      cout << "was here " << endl;
      data_scan[3*i]       = 0.0;
      data_scan[3*i + 1]   = 0.0;
      data_scan[3*i + 2]   = 0.0;

      colors_scan[3*i]     = 0;                     // r
      colors_scan[3*i+1]   = 0;                     // g
      colors_scan[3*i+2]   = 0;                     // b

      intensity_scan[i] = 0;
    }
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
      max_tmp = 32767;
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
      max_tmp = 32767;
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
          else if (intensity > cLevel2)
          {
            colors_scan[3*linesFile]     = 100;                         // r
            colors_scan[3*linesFile+1]   = 100;                     // g
            colors_scan[3*linesFile+2]   = 100;                         // b

            int abs_id = abs(id-512);
            colors_scan[3*abs_id ]     = 100;                         // r
            colors_scan[3*abs_id +1]   = 100;                     // g
            colors_scan[3*abs_id +2]   = 100;                         // b

          }
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
      max_tmp = 32767;
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
        filter(data_scan, colors_scan, intens);

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

int main(int argc, char* argv[])
{

/* Init */
  // Load parameters
  init(argc, argv);
  cout << "init finished" << endl;
  // Initialization of 3D viewer
  cloud_scan = new VtkCloud();
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
   boundryBox();

/* Create objects */
   initShowCloud();
   createShowCloud();
   // Create mirror plane
   // size = data_size_Scanner * 3_koordinates
   data_mirror   = new double[cloudsize_scan * 3];
   colors_mirror = new unsigned char[cloudsize_scan * 3];
   // Show sensor cube
   createCube(sensorsize);
   cout << "Cube created" << endl;
   // Show axes
   viewer3D->showAxes(true);
   cout << "Axes created" << endl;

/* Load and start vierwer */
   viewer3D->addCloud(cloud_show);
   viewer3D->addCloud(cloud_mirror);
   viewer3D->addCloud(cloud_sensor);
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

   delete cloud_scan;
   delete cloud_show;
   delete cloud_mirror;
   delete cloud_sensor;
   return 0;

}
