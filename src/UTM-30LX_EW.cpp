/**
* @file   ${UTM-30LX_EW.cpp}
* @author Rainer Koch
* @date   ${2014/01/09}
*
* Class for Hokuyo UTM-30LX_EW Laser scanner
*/

#define DEG2RAD M_PI/180.0

#include "Urg_driver.h"
#include "obcore/base/Logger.h"
#include "obcore/math/mathbase.h"
#include <iostream>
#include <stdio.h>
#include <csignal>
#include <cstdio>

#include <pthread.h>
#include <string.h>

#include "UTM-30LX_EW.h"

#define DIST_TH 10
#define INT_TH  10

using namespace laserDev;
using namespace qrk;


UTM30LXEW* _this;

pthread_t _thread;
bool _shutdown = false;
bool _run      = true;
pthread_mutex_t _mutex = PTHREAD_MUTEX_INITIALIZER;

void* taskCode(void*)
{
  while(_run)
  {
    _this->schedule();
    usleep(3000);
    //  pthread_yield();
  }

  _shutdown = true;
  //pthread_exit(NULL);

  return NULL;
}

UTM30LXEW::UTM30LXEW()
{
  _this = this;

  // Scanner parameters
  char* dev           = (char*)"192.168.0.10";
  long port_rate = 10940;  //115200 for Serial
  Lidar::connection_type_t dev_type = Lidar::Ethernet;
  Lidar::measurement_type_t meas_type = Lidar::Distance;
  int first_step = 135;
  int last_step = -135;
  int skip_step = -1;
  int measure_type = 0;
  int opendev = 0;


  std::vector<long> dist;
  std::vector<unsigned short> intens;

  /* Config with USB or Ethernet*/
     if(!(opendev = _laser.open(dev, port_rate, dev_type))){
       cout << "Can`t connect do device" << endl;
     }
     else
     {
       cout << "Connect to ip/dev: " << dev << " port/baudrate: " << port_rate << endl;

  /* Set scanner parameter */
       _laser.set_scanning_parameter(first_step, last_step, skip_step);
       cout << "Min step: " << _laser.min_step() << " Max step: " << _laser.max_step() << endl;
       cout << "Min dist.: " << _laser.min_distance() << " Max dist.: " << _laser.max_distance() << endl;
       cout << "Min echo size: " << _laser.max_echo_size() << " Max data size: " << _laser.max_data_size() << endl;

  /* Start scanning */
       _laser.laser_on();
       _laser.start_measurement(meas_type, -1, 0);


_nrOfRays = 1081;
/*TODO:
    if (_cfg.angleResolution == 2500) {
      _nrOfRays = 1081;
    }
    else if (_cfg.angleResolution == 5000) {
      _nrOfRays = 541;
    }
    else
    {
      LOGMSG(DBG_ERROR, "Unsupported resolution");
      exit(1);
    }*/

    _ranges       = new double[_nrOfRays];
    _intensities  = new double[_nrOfRays];
    _angles       = new double[_nrOfRays];
    _coords2D     = new double[2*_nrOfRays];
    _normals      = new double[2*_nrOfRays];
    _mask         = new bool  [_nrOfRays];



    const char* stat;
    do // wait for ready status
    {
      stat = _laser.status();
      sleep(1.0);
      cout << "waiting ... stat = " << stat << endl;
    }
    while (stat != 0);

    cout << "go" << endl;
//    _laser.scanContinous(1);

    pthread_create(&_thread, NULL, &taskCode, NULL);

  this->calculateAngles();
     }
}

UTM30LXEW::~UTM30LXEW()
{
  _run = false;
  //pthread_join( _thread, NULL);
  while(!_shutdown)
    usleep(100);

   cout << "stopping continuous mode" << endl;
   _laser.laser_off();

   cout << "disconnect" << endl;
   _laser.close();
}


double UTM30LXEW::getStartAngle()
{
  return _laser.index2deg(_laser.min_step());
}

double UTM30LXEW::getStopAngle()
{
  return _laser.index2deg(_laser.max_step());
}

unsigned int UTM30LXEW::getNumberOfRays()
{
  return _laser.max_data_size();
}

double UTM30LXEW::getAngularRes(void)
{
  return (-(_laser.index2deg(_laser.min_step()))+(_laser.index2deg(_laser.max_step())))/(_laser.max_data_size());
}

double* UTM30LXEW::getRanges()
{
  return _ranges;
}

double* UTM30LXEW::getCoords2D()
{
  return _coords2D;
}

void UTM30LXEW::schedule()
{
  pthread_mutex_lock(&_mutex);
  _laser.get_distance(dist);
  pthread_mutex_unlock(&_mutex);
}

bool UTM30LXEW::grab(void)
{
  if(_laser.is_stable())
  {
    pthread_mutex_lock(&_mutex);
    memcpy(&_distBuffer, &dist, sizeof(_distBuffer));
    pthread_mutex_unlock(&_mutex);

    this->calculateRanges();
    this->calculateIntensities();
    this->calculateCoords2D();
  }
  else
  {
    LOGMSG(DBG_ERROR, "Connection error");
  }

  return(true);
}

void UTM30LXEW::calculateRanges(void)
{
  _laser.get_distance(dist);
  cout << getNumberOfRays() << endl;
  for(unsigned int i=0; i<(unsigned int) _laser.max_data_size(); i++)
    _ranges[i] = dist[i] * 0.001;
}

void UTM30LXEW::calculateIntensities(void)
{
  //TODO:
  cout << " have to be done" << endl;
  //for(unsigned int i=0; i<(unsigned int)_data.rssi_len1; i++)
  //  _intensities[i] = _data.rssi1[i];
}

void UTM30LXEW::calculateAngles(void)
{
  double res = getAngularRes();
  for(unsigned int i=0 ; i<_nrOfRays ; i++)
    _angles[i] = getStartAngle() + ((double)i)*res;
}

void UTM30LXEW::calculateCoords2D(void)
{
  //TODO:
  for(int id=0; id <= getNumberOfRays(); id++)
  {
    // Build vtkCloud
    _coords2D[2*id]       = cos(_laser.index2rad(id)) * dist[id];  // x
    _coords2D[2*id + 1]   = sin(_laser.index2rad(id)) * dist[id];  // y
  }
}

