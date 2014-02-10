/**
* @file   ${hokuyoShow.cpp}
* @author Rainer Koch
* @date   ${2013/12/06}
*
* Stream Hokuyo to file, show Hokuyo measurements
* does use the standard urg driver
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <stdlib.h>
#include <math.h>
#include "Lidar.h"
#include "Urg_driver.h"
//#include "VisualObjects.h"

#include "hokuyoShow.h"
#include "obgraphic/Obvious3D.h"
#include <vector>

#define _USE_MATH_DEFINES

using namespace std;
using namespace qrk;
using namespace obvious;

    // 0 = Laser, 1 = File
   int sourcetype = 0;
/* Default parameters for laser*/
   char* dev           = (char*)"192.168.0.10";
   long port_rate = 10940;  //115200 for Serial
   Lidar::connection_type_t dev_type = Lidar::Ethernet;
   Lidar::measurement_type_t meas_type = Lidar::Distance;
   int first_step = 135;
   int last_step = -135;
   int skip_step = -1;
   int measure_type = 1;
   int opendev = 0;
   Urg_driver urg;

/* Default parameters for file*/
   char* file;
   unsigned int file_nr = 0;
   bool _save = true;

 /* Default parameters for viewer */
   std::vector<long> dist;
   std::vector<unsigned short> intens;
   VtkCloud* cloud_scan;
   VtkCloud* cloud_mirror;
   VtkCloud* cloud_sensor;

   int cloudsize;
   int raysPscan;
   int sensorsize = 5;  // only odd numbers
   double* data_scan;
   double* data_mirror;
   double* data_sensor;

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
  if(argc>1){
     int i = atoi(argv[1]);
     if(i==0) dev_type = Lidar::Serial;
     else if(i==1) dev_type = Lidar::Ethernet;
     else if(i==2)
     {
       file = argv[2];
       sourcetype = 1;
       cout << "Load file" << argv[2] << endl;
     }
     else cout << "No device number! 0=Ethernet, 1=Serial" << endl;
   }

   if(argc>2)
   {
      if(argc!=5)
      {
         cout << "usage: " << argv[0] << " [type (0=Serial 1=Ethernet 2=File)] [ip/dev/filename] [port/baudrate/0] [measurements(0=data, 1=data+intensity, 2=multi_echo_data, 3=multi_echo_data+intensities)" << endl;
      }
      dev           = argv[2];
      port_rate     = atoi(argv[3]);
      measure_type  = atoi(argv[4]);

   }

   switch(measure_type)
    {
    case 0:
      meas_type =  Lidar::Distance;
      cout << "Distance" << endl;
      break;
    case 1:
      meas_type =  Lidar::Distance_intensity;
      cout << "Distance + Intensity" << endl;
      break;
    case 2:
      meas_type =  Lidar::Multiecho;
      cout << "Multiecho" << endl;
      break;
    case 3:
      meas_type =  Lidar::Multiecho_intensity;
      cout << "Multiecho + Intensity" << endl;
      break;
    default:
      meas_type =  Lidar::Distance;
      cout << "Default: Distance" << endl;
    }
}
void save_raw_dist(std::vector<long>& distance)
{
    // Save as raw file
    ofstream file;
    char filename[64];
    static int cnt = 0;

    double degree = 0.0;
    sprintf(filename, "scan_raw_%03d.txt", cnt);
    file.open (filename);
    file << "angle, distance" << endl;

    for(int i=0; i<= raysPscan; i++)
    {
      degree = urg.index2rad(i);
      file << degree << " " << distance[i] << endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}
void save_xy_dist(std::vector<long>& distance)
{
    // Save as xy file
    ofstream file;
    char filename[64];
    static int cnt = 0;
    double x = 0.0;
    double y = 0.0;
    sprintf(filename, "scan_xy_%03d.txt", cnt);
    file.open (filename);
   file << "x, y" << endl;

    for(int i=0; i<= raysPscan; i++)
    {
      x = sin(urg.index2rad(i)) * dist[i];
      y = cos(urg.index2rad(i)) * dist[i];
      file << "" << x << " " << y << endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}
void save_raw_multidist(std::vector<long>& distance)
{
    // Save as raw file
    ofstream file;
    char filename[64];
    static int cnt = 0;
    unsigned int id = 0;
    double degree = 0.0;
    sprintf(filename, "scan_raw_m_%03d.txt", cnt);
    file.open (filename);
    file << "angle, distance1, distance2, distance3" << endl;
    for(int i=0; i<= raysPscan; i++)
    {
      degree = urg.index2rad(i);
      file << degree << " " << distance[3*i] << " " << distance[3*i+1] << " " << distance[3*i+2] << endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}
void save_xy_multidist(std::vector<long>& distance)
{
    // Save as xy file
    ofstream file;
    char filename[64];
    static int cnt = 0;
    double x = 0.0;
    double y = 0.0;
    double degree = 0.0;
    sprintf(filename, "scan_xy_m_%03d.txt", cnt);
    file.open (filename);
    file << "angle, x, y " << endl;

    for(int i=0; i<= raysPscan; i++)
    {
        degree = urg.index2rad(i);
        x = sin(urg.index2rad(i)) * dist[3*i];
        y = cos(urg.index2rad(i)) * dist[3*i];
        file << degree << " " << x << " " << y <<  endl;
        x = sin(urg.index2rad(i)) * dist[3*i+1];
        y = cos(urg.index2rad(i)) * dist[3*i+1];
        file << degree << " " << x << " " << y << endl;
        x = sin(urg.index2rad(i)) * dist[3*i+2];
        y = cos(urg.index2rad(i)) * dist[3*i+2];
        file << degree << " " << x << " " << y << endl;
    }


//    while(id <= 3*urg.max_data_size())
//    {
//      degree = urg.index2rad(id);
//      x = cos(urg.index2rad(id)) * dist[id];
//      y = sin(urg.index2rad(id)) * dist[id];
//      file << degree << " " << x << " " << y << endl;
//      x = cos(urg.index2rad(id)) * dist[id+1];
//      y = sin(urg.index2rad(id)) * dist[id+1];
//      file << degree << " " << x << " " << y << endl;
//      x = cos(urg.index2rad(id)) * dist[id+2];
//      y = sin(urg.index2rad(id)) * dist[id+2];
//      file << degree << " " << x << " " << y << endl;
//
//      id = id + 3;
//    }
    file.close();
    cout << "saved: " << filename << endl;

    cnt++;
}
void save_raw_intens(std::vector<long>& distance, std::vector<unsigned short>& intensity)
{
    // Save as raw file
    ofstream file;
    char filename[64];
    static int cnt = 0;

    double degree = 0.0;
    sprintf(filename, "scan_raw_i_%03d.txt", cnt);
    file.open (filename);
    file << "angle, distance, intensity" << endl;
    for(int i=0; i<= raysPscan; i++)
    {
      degree = urg.index2rad(i);
      file << degree << " " << distance[i] << " " << intensity[i] << endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}
void save_xy_intens(std::vector<long>& distance, std::vector<unsigned short>& intensity)
{
    // Save as xy file
    ofstream file;
    char filename[64];
    static int cnt = 0;
    double x = 0.0;
    double y = 0.0;
    sprintf(filename, "scan_xy_i_%03d.txt", cnt);
    file.open (filename);
    file << "x, y, intensity" << endl;
    for(int i=0; i<= raysPscan; i++)
    {
      x = sin(urg.index2rad(i)) * dist[i];
      y = cos(urg.index2rad(i)) * dist[i];
      file << x << " " << y << " " << intensity[i] << endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}
void save_raw_mulitintens(std::vector<long>& distance, std::vector<unsigned short>& intensity)
{

    // Save as raw file
    ofstream file;
    char filename[64];
    static int cnt = 0;
    double degree = 0.0;
    sprintf(filename, "scan_raw_mi_%03d.txt", cnt);
    file.open (filename);
    file << "angle, distance0, distance1, distance2, intensity0, intensity1, intensity2" << endl;

    for(int i=0; i<= raysPscan; i++)
    {
       degree = urg.index2rad(i);
       file << degree << " " << distance[i] << " " << distance[3*i+1] << " " << distance[3*i+2] << " " << intensity[3*i] << " " << intensity[3*i+1] << " " << intensity[3*i+2] << endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}
void save_xy_mulitintens(std::vector<long>& distance, std::vector<unsigned short>& intensity)
{

    // Save as xy file
    ofstream file;
    char filename[64];
    static int cnt = 0;
    double x = 0.0;
    double y = 0.0;
    double degree = 0.0;
    sprintf(filename, "scan_xy_mi_%03d.txt", cnt);
    file.open (filename);
    file << "angle, x, y, intensity " << endl;
    for(int i=0; i<= raysPscan; i++)
    {
      degree = urg.index2rad(i);
      x = sin(urg.index2rad(i)) * dist[3*i];
      y = cos(urg.index2rad(i)) * dist[3*i];
      file << degree << " " << x << " " << y << " " << intensity[i]<< endl;
      x = sin(urg.index2rad(i)) * dist[3*i+1];
      y = cos(urg.index2rad(i)) * dist[3*i+1];
      file << degree << " " << x << " " << y << " " << intensity[i+1]<< endl;
      x = sin(urg.index2rad(i)) * dist[3*i+2];
      y = cos(urg.index2rad(i)) * dist[3*i+2];
      file << degree << " " << x << " " << y << " " << intensity[i+2]<< endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}

void createCupe(int cubesize)
{
  cout << "be here" << endl;
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
  cout << "be here done" << endl;

}
void startLaser()
{
  /* Config with USB or Ethernet*/
  if(!(opendev = urg.open(dev, port_rate, dev_type)))
  {
   cout << "Can`t connect do device" << endl;
  }
  else
  {
   cout << "Connect to ip/dev: " << dev << " port/baudrate: " << port_rate << endl;

   /* Set scanner parameter */
   urg.set_scanning_parameter(first_step, last_step, skip_step);
   urg.set_measurement_type(meas_type);

   cout << "Min step: " << urg.min_step() << " Max step: " << urg.max_step() << endl;
   cout << "Min dist.: " << urg.min_distance() << " Max dist.: " << urg.max_distance() << endl;
   cout << "Min echo size: " << urg.max_echo_size() << " Max data size: " << urg.max_data_size() << endl;

   /* Start scanning */
   urg.laser_on();
   urg.start_measurement(meas_type, -1, 0);
  }
}
void stopLaser()
{
  urg.stop_measurement();
  urg.laser_off();
  urg.close();
}

int load_xy_file(char* filename)
{
  // Read X,Y file
      float x, y, degree, intensity;
      string leseString;
      ifstream file(filename);
      int echoCount = 0;
      max_tmp = 32767;
      unsigned int linesFile;
      unsigned int id;

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
              data_scan[3*linesFile + 2]   = 0;                                  // z

              intens[linesFile] = intensity;

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
              //cout << "data_scan id: " << linesFile << " x: " << data_scan[3*linesFile] << "y: " << data_scan[3*linesFile+1] << endl;
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
              intens[linesFile] = intensity;

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
void cbSave()
{
  std::vector<long> distance;
  std::vector<unsigned short> intensity;
  /* Save file */
       switch(measure_type)
       {
       case 0:
         cout << "distance save" << endl;
         urg.get_distance(distance);
         save_raw_dist(distance);
         save_xy_dist(distance);
         cout << "distance finish " << endl;
         break;
       case 1:
         cout << "distance + intensity save" << endl;
         urg.get_distance_intensity(distance, intensity);
         save_raw_intens(distance, intensity);
         save_xy_intens(distance, intensity);
         cout << "distance + intensity finish " << endl;
         break;
       case 2:
         cout << "multi distance save" << endl;
         urg.get_multiecho(distance);
         save_raw_multidist(distance);
         save_xy_multidist(distance);
         cout << "multi distance finish " << endl;
         break;
       case 3:
         cout << "multi distance + intensity save" << endl;
         urg.get_multiecho_intensity(distance, intensity);
         save_raw_mulitintens(distance, intensity);
         save_xy_mulitintens(distance, intensity);
         cout << "multi distance + intensity finish" << endl;
         break;
       default:
         cout << "Default: Multiecho + Intensity" << endl;
       }
}
void cbMeasuretype()
{
  if(measure_type < 3)
  {
    measure_type++;
  }
  else
  {
    measure_type = 0;
  }

   // delete old data
  for(int i=0; i <= cloudsize; i++)
           {
             // Build vtkCloud
             data_scan[3*i]       = 0;  // x
             data_scan[3*i + 1]   = 0;  // y
             data_scan[3*i + 2]   = 0;                                  // z

             colors_scan[3*i]     = 0;                                // r
             colors_scan[3*i+1]   = 0;                                // g
             colors_scan[3*i+2]   = 0;                                // b
           }

  switch(measure_type)
   {
   case 0:
     stopLaser();
     meas_type =  Lidar::Distance;
     startLaser();
     cout << "Distance" << endl;
     break;
   case 1:
     stopLaser();
     meas_type =  Lidar::Distance_intensity;
     startLaser();
     cout << "Distance + Intensity" << endl;
     break;
   case 2:
     stopLaser();
     meas_type =  Lidar::Multiecho;
     startLaser();
     cout << "Multiecho" << endl;
     break;
   case 3:
     //TODO: Not finished
     stopLaser();
     meas_type =  Lidar::Multiecho_intensity;
     startLaser();
     cout << "Multiecho + Intensity" << endl;
     break;
   }
  cout << "Cloudsize: " << cloudsize << endl;
  test = 0;

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
    if(!_pause)
    {
      int i= 0;
/* Show data */
      if(sourcetype == 0)
      {
        switch(measure_type)
         {
         case 0:
           urg.get_distance(dist);

           for(int i=0; i <= raysPscan; i++)
           {
             // Build vtkCloud
             data_scan[3*i]       = sin(urg.index2rad(i)) * dist[i];  // x
             data_scan[3*i + 1]   = cos(urg.index2rad(i)) * dist[i];  // y
             data_scan[3*i + 2]   = 0;                                  // z

             colors_scan[3*i]     = 200;                                // r
             colors_scan[3*i+1]   = 200;                                // g
             colors_scan[3*i+2]   = 200;                                // b
           }
           break;
         case 1:
           urg.get_distance_intensity(dist, intens);

           for(int i=0; i <= raysPscan; i++)
           {
  //TODO: change colors in better way
             // Change color on the intensity value
             // intensColor = intens[i]/max_tmp*256;

             // Change color by limits
             if(intens[i]<= 2700)
             {
               colors_scan[3*i]     = 100;                         // r
               colors_scan[3*i+1]   = 100;                     // g
               colors_scan[3*i+2]   = 100;                         // b
             }
             else
             {
               colors_scan[3*i]     = 255;                         // r
               colors_scan[3*i+1]   = 100;                     // g
               colors_scan[3*i+2]   = 0;                         // b
             }

             // Build vtkCloud
             data_scan[3*i]       = sin(urg.index2rad(i)) * dist[i];  // x
             data_scan[3*i + 1]   = cos(urg.index2rad(i)) * dist[i];  // y
             data_scan[3*i + 2]   = intens[i]/10;                     // z
           }
           break;
         case 2:
           urg.get_multiecho(dist);
           for(int i=0; i <= raysPscan; i++)
           {
             // Build vtkCloud
               // First echo
               data_scan[9*i]       = sin(urg.index2rad(i)) * dist[3*i];    // x
               data_scan[9*i + 1]   = cos(urg.index2rad(i)) * dist[3*i];    // y
               data_scan[9*i + 2]   = 0;                                    // z

               colors_scan[9*i]     = 200;                                  // r
               colors_scan[9*i+1]   = 200;                                  // g
               colors_scan[9*i+2]   = 200;                                  // b

               // Second echo
               data_scan[9*i + 3]   = sin(urg.index2rad(i)) * dist[3*i+1];  // x
               data_scan[9*i + 4]   = cos(urg.index2rad(i)) * dist[3*i+1];  // y
               data_scan[9*i + 5]   = 25;                                   // z

               colors_scan[9*i+3]   = 200;                                  // r
               colors_scan[9*i+4]   = 0;                                    // g
               colors_scan[9*i+5]   = 0;                                    // b

               // Third echo
               data_scan[9*i + 6]   = sin(urg.index2rad(i)) * dist[3*i+2];  // x
               data_scan[9*i + 7]   = cos(urg.index2rad(i)) * dist[3*i+2];  // y
               data_scan[9*i + 8]   = 50;                                   // z

               colors_scan[9*i+6]   = 0;                                    // r
               colors_scan[9*i+7]   = 200;                                    // g
               colors_scan[9*i+8]   = 0;                                  // b
           }
           break;
         case 3:
            urg.get_multiecho_intensity(dist, intens);

            for(int i=0; i <= raysPscan; i++)
            {

               if (intens[i] < min_tmp)
                 min_tmp = intens[i];
               if (intens[i] > max_tmp)
                 max_tmp = intens[i];

//                if(test1==0)
//                cout << intens[i] << "  " ;
            }
//            if(test1==0)
//              cout << endl;
//            cout << "Min: " << min_tmp << " Max: " << max_tmp << endl;
//            test1 = 1;

            for(int i=0; i <= raysPscan; i++)
            {
             //TODO: change colors in better way
             intensColor = intens[i]/max_tmp*256;

             // Build vtkCloud
                // First echo
                data_scan[9*i]       = sin(urg.index2rad(i)) * dist[3*i];    // x
                data_scan[9*i + 1]   = cos(urg.index2rad(i)) * dist[3*i];    // y
                data_scan[9*i + 2]   = intens[3*i]/10;                                    // z

                colors_scan[9*i]     = intensColor;                         // r
                colors_scan[9*i+1]   = 100;                     // g
                colors_scan[9*i+2]   = 256-intensColor;                         // b

                // Second echo
                data_scan[9*i + 3]   = sin(urg.index2rad(i)) * dist[3*i+1];  // x
                data_scan[9*i + 4]   = cos(urg.index2rad(i)) * dist[3*i+1];  // y
                data_scan[9*i + 5]   = 25;                                   // z

                colors_scan[9*i+3]   = 256-intensColor;                         // r
                colors_scan[9*i+4]   = intensColor;                     // g
                colors_scan[9*i+5]   = 100;                         // b

                // Third echo
                data_scan[9*i + 6]   = sin(urg.index2rad(i)) * dist[3*i+2];  // x
                data_scan[9*i + 7]   = cos(urg.index2rad(i)) * dist[3*i+2];  // y
                data_scan[9*i + 8]   = 50;                                   // z

                colors_scan[9*i+6]   = 100;                         // r
                colors_scan[9*i+7]   = 256-intensColor;                     // g
                colors_scan[9*i+8]   = intensColor;                         // b
            }
         break;
         }
      }

/* Filter */

      float distance_refl = 0.0;
      float angle_refl = 0.0;
      float angle_refl_g = 0.0;

      float distance_orig = 0.0;
      float angle_orig = 0.0;
      float angle_orig_g = 0.0;

      unsigned int id_orig = 0;

      if(_filterSwitch && measure_type == 1)
      {
        for(int i=0; i<=raysPscan; i++)
        {
          // threshold for intensity
          if(intens[i] > 2700)
          {
            // recalculate original point id
            id_orig = abs(raysPscan-i);
 //           cout << "id_refl: " << i << " id_orig: " << id_orig << endl;

            // calculate vector
            //TODO: later expand with z => sqrt(x2+y2+z2)
            distance_refl = sqrt(((data_scan[3*i])*(data_scan[3*i])) + (data_scan[3*i+1]*data_scan[3*i+1]));
            angle_refl_g = asin(data_scan[3*i]/distance_refl*1.0)/M_PI*360;
            angle_refl = asin(data_scan[3*i]/distance_refl*1.0);

            distance_orig = sqrt(((data_scan[3*id_orig])*(data_scan[3*id_orig])) + (data_scan[3*id_orig+1]*data_scan[3*id_orig+1]));
            angle_orig_g = abs(360.0-angle_refl);
            angle_refl = M_PI-asin(data_scan[3*i]/distance_refl*1.0);


            // mark reflected point
            colors_scan[3*i]     = 100;                         // r
            colors_scan[3*i +1]   = 0;                     // g
            colors_scan[3*i +2]   = 0;                         // b

            // mark original point
            colors_scan[3*id_orig]     = 100;                         // r
            colors_scan[3*id_orig +1]   = 0;                     // g
            colors_scan[3*id_orig +2]   = 100;                         // b

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
              cout << " x_refl: " << data_scan[3*i] << " y_ref: " << data_scan[3*i+1] << endl;
              cout << " x_orig: " << data_scan[3*id_orig] << " y_orig: " << data_scan[3*id_orig+1] << endl;
              cout << " x_mirr: " << data_mirror[3*i] << " y_mirr: " << data_mirror[3*i+1] << endl << endl << endl;
            }
          }
        }

        test = 1;

        /* Set recalculated mirro plane */
        cloud_mirror->setCoords(data_mirror, cloudsize, 3);
        cloud_mirror->setColors(colors_mirror, cloudsize, 3);

      }

/* Set new cloud*/
        cloud_scan->setCoords(data_scan, cloudsize, 3);
        cloud_scan->setColors(colors_scan, cloudsize, 3);

      viewer3D->update();
     }
    }

private:

};




int main(int argc, char* argv[])
{

/* Init */
  // Load parameters
  init(argc, argv);
  cout << "init finished" << endl;

/* Initialization of scanner */
  if(sourcetype==0)
  {
    startLaser();
    cout << "Laser started" << endl;
  }

/* Initialization of 3D viewer */
   cloud_scan = new VtkCloud();
   cloud_mirror = new VtkCloud();

   /* Show sensor */
   cloud_sensor = new VtkCloud();
   createCupe(sensorsize);
   /* Show axes */
   viewer3D->showAxes(true);

   // size = data_size_Scanner * 3_echos * 3_koordinates
   data_mirror   = new double[cloudsize * 3 * 3];
   colors_mirror = new unsigned char[cloudsize * 3 * 3];

   if(sourcetype == 0)
   {
     cloudsize = 3*urg.max_data_size();
     raysPscan = urg.max_data_size();
     data_scan   = new double[cloudsize * 3 * 3];
     colors_scan = new unsigned char[cloudsize * 3 * 3];
   }
   else if(sourcetype == 1)
   {
     cloudsize = load_xy_file(file);
   }

   // define update timer
   vtkSmartPointer<vtkTimerCallback> cb =  vtkSmartPointer<vtkTimerCallback>::New();
   vtkSmartPointer<vtkRenderWindowInteractor> interactor = viewer3D->getWindowInteractor();
   interactor->AddObserver(vtkCommand::TimerEvent, cb);
   interactor->CreateRepeatingTimer(30);
   // define keyboard callbacks
   viewer3D->registerFlipVariable("p", &_pause);
   viewer3D->registerKeyboardCallback("y", cbSave);
   viewer3D->registerKeyboardCallback("m", cbMeasuretype);
   viewer3D->registerFlipVariable("f", &_filterSwitch);

   cout << "VTK vierwer init" << endl;

   viewer3D->addCloud(cloud_scan);
   viewer3D->addCloud(cloud_mirror);
   viewer3D->addCloud(cloud_sensor);
   viewer3D->startRendering();


/* Shutdown */
   stopLaser();
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
