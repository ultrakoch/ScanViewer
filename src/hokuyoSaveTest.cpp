/**
* @file   ${hokuyoSaveTest.cpp}
* @author Rainer Koch
* @date   ${2013/12/06}
*
* Stream Hokuyo to file
* does use the standard urg driver
*/

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <stdlib.h>
#include <math.h>
#include "Lidar.h"
#include "Urg_driver.h"

#include "hokuyoSaveTest.h"
#include "obgraphic/Obvious2D.h"

using namespace std;
using namespace qrk;
using namespace obvious;


/* Default parameters */
   char* dev           = (char*)"192.168.0.10";
   long port_rate = 10940;  //115200 for Serial
   Lidar::connection_type_t dev_type = Lidar::Ethernet;
   Lidar::measurement_type_t meas_type = Lidar::Distance;
   int first_step = 135;
   int last_step = -135;
   int skip_step = -1;
   int measure_type = 0;
   int opendev = 0;
   Urg_driver urg;
   char* path;
   int file_nr = 0;
   bool _save = true;

   std::vector<long> dist;
   std::vector<unsigned short> intens;


void init(int argc, char* argv[])
{
  if(argc>1){
     int i = atoi(argv[1]);
     if(i==0) dev_type = Lidar::Serial;
     else if(i==1) dev_type = Lidar::Ethernet;
     else cout << "No device number! 0=Ethernet, 1=Serial" << endl;
   }

   if(argc>2)
   {
      if(argc!=5)
      {
         cout << "usage: " << argv[0] << " [type (0=Serial 1=Ethernet)] [ip/dev] [port/baudrate] [measurements(0=data, 1=data+intensity, 2=multi_echo_data, 3=multi_echo_data+intensities)" << endl;
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
      cout << "Default: Multiecho + Intensity" << endl;
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

    for(int id=0; id<= urg.max_data_size(); id++)
    {
      degree = urg.index2rad(id);
      file << degree << " " << distance[id] << endl;
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

    for(int id=0; id<= urg.max_data_size(); id++)
    {
      x = cos(urg.index2rad(id)) * dist[id];
      y = sin(urg.index2rad(id)) * dist[id];
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
//    while(id <= 3*urg.max_data_size())
//    {
//      degree = urg.index2rad(id);
//      file << degree << " " << distance[id] << " " << distance[id+1] << " " << distance[id+2] << endl;
//      id = id + 3;
//    }
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
    unsigned int id = 0;
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
    for(int id=0; id<= urg.max_data_size(); id++)
    {
      degree = urg.index2rad(id);
      file << degree << " " << distance[id] << " " << intensity[id] << endl;
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
    for(int id=0; id<= urg.max_data_size(); id++)
    {
      x = cos(urg.index2rad(id)) * dist[id];
      y = sin(urg.index2rad(id)) * dist[id];
      file << x << " " << y << " " << intensity[id] << endl;
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
    unsigned int id = 0;
    double degree = 0.0;
    sprintf(filename, "scan_raw_mi_%03d.txt", cnt);
    file.open (filename);
    file << "angle, distance0, distance1, distance2, intensity0, intensity1, intensity2" << endl;

    while(id <= 3*urg.max_data_size())
     {
       degree = urg.index2rad(id);
       file << degree << " " << distance[id] << " " << distance[id+1] << " " << distance[id+2] << " " << intensity[id] << " " << intensity[id+1] << " " << intensity[id+2] << endl;
       id = id + 3;
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
    unsigned int id = 0;
    double x = 0.0;
    double y = 0.0;
    double degree = 0.0;
    sprintf(filename, "scan_xy_mi_%03d.txt", cnt);
    file.open (filename);
    file << "angle, x, y, intensity " << endl;
    while(id <= 3*urg.max_data_size())
    {
      degree = urg.index2rad(id);
      x = cos(urg.index2rad(id)) * dist[id];
      y = sin(urg.index2rad(id)) * dist[id];
      file << degree << " " << x << " " << y << " " << intensity[id]<< endl;
      x = cos(urg.index2rad(id)) * dist[id+1];
      y = sin(urg.index2rad(id)) * dist[id+1];
      file << degree << " " << x << " " << y << " " << intensity[id+1]<< endl;
      x = cos(urg.index2rad(id)) * dist[id+2];
      y = sin(urg.index2rad(id)) * dist[id+2];
      file << degree << " " << x << " " << y << " " << intensity[id+2]<< endl;
      id = id + 3;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
}

int main(int argc, char* argv[])
{

/* Init */
  // Load parameters
  init(argc, argv);

/* Config with USB or Ethernet*/
   if(!(opendev = urg.open(dev, port_rate, dev_type))){
     cout << "Can`t connect do device" << endl;
   }
   else
   {
     cout << "Connect to ip/dev: " << dev << " port/baudrate: " << port_rate << endl;

/* Set scanner parameter */
     urg.set_scanning_parameter(first_step, last_step, skip_step);
          cout << "Min step: " << urg.min_step() << " Max step: " << urg.max_step() << endl;
          cout << "Min dist.: " << urg.min_distance() << " Max dist.: " << urg.max_distance() << endl;
          cout << "Min echo size: " << urg.max_echo_size() << " Max data size: " << urg.max_data_size() << endl;

/* Start scanning */
     urg.laser_on();
     urg.start_measurement(meas_type, -1, 0);

/* Save file */
     switch(measure_type)
     {
     case 0:
       cout << "distance save" << endl;
       urg.get_distance(dist);
       save_raw_dist(dist);
       save_xy_dist(dist);
       cout << "distance finish " << endl;
       break;
     case 1:
       cout << "distance + intensity save" << endl;
       urg.get_distance_intensity(dist, intens);
       save_raw_intens(dist, intens);
       save_xy_intens(dist, intens);
       cout << "distance + intensity finish " << endl;
       break;
     case 2:
       cout << "multi distance save" << endl;
       urg.get_multiecho(dist);
       save_raw_multidist(dist);
       save_xy_multidist(dist);
       cout << "multi distance finish " << endl;
       break;
     case 3:
       cout << "multi distance + intensity save" << endl;
       urg.get_multiecho_intensity(dist, intens);
       save_raw_mulitintens(dist, intens);
       save_xy_mulitintens(dist, intens);
       cout << "multi distance + intensity finish" << endl;
       break;
     default:
       meas_type =  Lidar::Distance;
       cout << "Default: Multiecho + Intensity" << endl;
     }

/* Shutdown */
   urg.laser_off();
   urg.close();
   return 0;
   }
}
