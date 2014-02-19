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
/* Default parameters for laser*/
//   char* dev           = (char*)"192.168.0.10";
//   long port_rate = 10940;  //115200 for Serial
//   Lidar::connection_type_t dev_type = Lidar::Ethernet;
//   Lidar::measurement_type_t meas_type = Lidar::Distance;
//   int first_step = 135;
//   int last_step = -135;
//   int skip_step = -1;
//   int measure_type = 1;
//   int opendev = 0;
//   Urg_driver urg;

/* Default parameters for file*/
   char* file;
   string filetype;
   int measure_type;
   unsigned int file_nr = 0;
   bool _save = true;

 /* Default parameters for viewer */
   std::vector<long> dist;
   float* intens;
   VtkCloud* cloud_scan;
   VtkCloud* cloud_mirror;
   VtkCloud* cloud_sensor;

   int cloudsize;
   int raysPscan;
   int sensorsize = 5;  // only odd numbers
   double* data_scan;
   double* data_mirror;
   double* data_sensor;
   double* intensity_scan;

   double* sensorPos;
   unsigned char* colors_scan;
   unsigned char* colors_mirror;
   unsigned char* colors_sensor;
   unsigned char intensColor;

   Obvious3D* viewer3D = new Obvious3D();
   bool _pause = false;
   bool _filterSwitch = false;
   int id = 0;            // hokuyo measurement id

   // Testvariables
   bool test = 0;
   bool test1 = 0;
   unsigned int max_tmp = 0;
   unsigned int min_tmp = 65535;
   unsigned int tmp = 0;

void init(int argc, char* argv[])
{
  if(argc<3)
  {
     cout << "usage: " << argv[0] << "[filename] [txt_xy, rxp ] [measurements(0=data, 1=data+intensity, 2=multi_echo_data, 3=multi_echo_data+intensities)]" << endl;
  }
  else
  {
   file = argv[1];
   filetype = argv[2];
   measure_type = atoi(argv[3]);

   cout << "Load file" << file << endl;


   switch(measure_type)
   {
    case 0:
      cout << "Distance" << endl;
      break;
    case 1:
      cout << "Distance + Intensity" << endl;
      break;
    case 2:
      cout << "Multiecho" << endl;
      break;
    case 3:
      cout << "Multiecho + Intensity" << endl;
      break;
   }
  }
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


    for(int i=0; i<=raysPscan; i++)
    {
      // threshold for intensity
      if(intensity[i] > 2700)
      {
        // recalculate original point id
        id_orig = abs(raysPscan-i);
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
    cloud_mirror->setCoords(data_mirror, cloudsize, 3);
    cloud_mirror->setColors(colors_mirror, cloudsize, 3);

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

        switch(measure_type)
        {
        case 0:
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
            //cout << "lines: " << linesFile << endl;
            linesFile = 0;
            data_scan   = new double[id * 3]; // linesFile - header
            colors_scan = new unsigned char[id * 3];
            std::ifstream file(filename);
            // print file header and remove it
            getline(file, leseString);
          //  cout << "Fileheader:" << leseString << endl;
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
              //cout << "data id: " << linesFile << " x: " << data[3*linesFile] << "y: " << data[3*linesFile+1] << endl;
              linesFile++;
             }
          }
          break;
        case 1:
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


              if(intensity<= 3000)
              {
                colors_scan[3*linesFile]     = 100;                         // r
                colors_scan[3*linesFile+1]   = 0;                     // g
                colors_scan[3*linesFile+2]   = 0;                         // b
              }
              else if ((intensity > 3000) && (intensity<= 4000))
              {
                colors_scan[3*linesFile]     = 0;                         // r
                colors_scan[3*linesFile+1]   = 100;                     // g
                colors_scan[3*linesFile+2]   = 0;                         // b
              }
              else if (intensity > 4000)
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
          break;
        case 2:
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
          break;
        case 3:
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
          break;
        }
        return id;
}

int load_rxp_file(char* filename)
{
  // Read rxp file
      float x, y, z, degree, intensity, amplitude, deviation;
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
         cout << "lines: " << linesFile << endl;
         linesFile = 0;
         data_scan   = new double[id * 3]; // linesFile - header
         intensity_scan   = new double[id]; // linesFile - header
         colors_scan = new unsigned char[id * 3];
         std::ifstream file(filename);
         // print file header and remove it
         getline(file, leseString);
         cout << "Fileheader:" << leseString << endl;
         cout << "Fileheader:" << leseString << endl;

         for(leseString; getline(file, leseString);)
         {
           istringstream(leseString,ios_base::in) >> x >> y >> z >> intensity >> amplitude >> deviation;

           // Build vtkCloud
           data_scan[3*linesFile]       = x;
           data_scan[3*linesFile + 1]   = y;
           data_scan[3*linesFile + 2]   = z;

           intensity_scan[linesFile] = intensity;

// TODO: Amplitude and deviation not used yet

           if(intensity<= 3000)
           {
             colors_scan[3*linesFile]     = 100;                         // r
             colors_scan[3*linesFile+1]   = 0;                     // g
             colors_scan[3*linesFile+2]   = 0;                         // b
           }
           else if ((intensity > 3000) && (intensity<= 4000))
           {
             colors_scan[3*linesFile]     = 0;                         // r
             colors_scan[3*linesFile+1]   = 100;                     // g
             colors_scan[3*linesFile+2]   = 0;                         // b
           }
           else if (intensity > 4000)
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
      if(_filterSwitch && measure_type == 1)
        filter(data_scan, colors_scan, intens);

/* Set new cloud*/
        cloud_scan->setCoords(data_scan, cloudsize, 3);
        cloud_scan->setColors(colors_scan, cloudsize, 3);

      viewer3D->update();

    }

private:

};




int main(int argc, char* argv[])
{

/* Init */
  // Load parameters
  init(argc, argv);
  cout << "init finished" << endl;


/* Initialization of 3D viewer */
   cloud_scan = new VtkCloud();
   cloud_mirror = new VtkCloud();

   /* Show sensor */
   cloud_sensor = new VtkCloud();
   createCube(sensorsize);
   cout << "Cube created" << endl;

   /* Show axes */
   viewer3D->showAxes(true);
   cout << "Axes created" << endl;

   // size = data_size_Scanner * 3_echos * 3_koordinates
   data_mirror   = new double[cloudsize * 3 * 3];
   colors_mirror = new unsigned char[cloudsize * 3 * 3];

   if(filetype.compare("txt_xy"))
     cloudsize = load_xy_file(file);
   else if(filetype.compare("rxp"))
     cloudsize = load_rxp_file(file);
   else
     cout << "Can not load filetype " << filetype << endl;

   if (cloudsize!=0)
   cout << "File " << file << " loaded" << endl;


   // define update timer
   vtkSmartPointer<vtkTimerCallback> cb =  vtkSmartPointer<vtkTimerCallback>::New();
   vtkSmartPointer<vtkRenderWindowInteractor> interactor = viewer3D->getWindowInteractor();
   interactor->AddObserver(vtkCommand::TimerEvent, cb);
   interactor->CreateRepeatingTimer(30);
   // define keyboard callbacks
   viewer3D->registerFlipVariable("f", &_filterSwitch);

   cout << "VTK vierwer init" << endl;

   viewer3D->addCloud(cloud_scan);
   viewer3D->addCloud(cloud_mirror);
   viewer3D->addCloud(cloud_sensor);
   viewer3D->startRendering();


/* Shutdown */
   delete viewer3D;

   delete [] data_scan;
   delete [] data_mirror;
   delete [] data_sensor;

   delete [] colors_scan;
   delete [] colors_mirror;
   delete [] colors_sensor;

   delete cloud_scan;
   delete cloud_mirror;
   delete cloud_sensor;
   return 0;

}
