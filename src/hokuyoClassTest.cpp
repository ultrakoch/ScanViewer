/**
* @file   ${hokuyoClassTest.cpp}
* @author Rainer Koch
* @date   ${2013/12/06}
*
* Stream Hokuyo to file
* use the own Laser class
*/

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <stdlib.h>
#include <math.h>
#include "UTM-30LX_EW.h"

#include "hokuyoClassTest.h"
#include "obgraphic/Obvious3D.h"
#include <vector>

using namespace laserDev;
using namespace std;
using namespace obvious;

UTM30LXEW* hokuyo;

/* Default parameters */

   char* path;
   int file_nr = 0;
   bool _save = true;

   double* dist;
   std::vector<unsigned short> intens;


   VtkCloud* cloud;
   double* data;
   unsigned char* colors;
   unsigned char intensColor;
   int cloudsize;
   int measure_type = 0;
   Obvious3D* viewer3D = new Obvious3D();
   bool              _pause       = false;

void init(int argc, char* argv[])
{

}
void save_raw_dist(double*  distance, int rays)
{
    // Save as raw file
    ofstream file;
    char filename[64];
    static int cnt = 0;

    double degree = 0.0;
    sprintf(filename, "scan_raw_%03d.txt", cnt);
    file.open (filename);
    file << "angle, distance" << endl;

    for(int id=0; id<= rays; id++)
    {
 //TODO:     degree = _laser.index2rad(id);//urg.index2rad(id);
      file << degree << " " << distance[id] << endl;
    }
    file.close();
    cout << "saved: " << filename << endl;
    cnt++;
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
      switch(measure_type)
       {
       case 0:
        hokuyo->grab();
        data = hokuyo->getCoords2D();
         for(int id=0; id < 1081; id++)
         {
           // Build vtkCloud
           cloud[3*id]       = data[2*id];  // x
           cloud[3*id + 1]   = data[2*id+1];  // y
           cloud[3*id + 2]   = 0;                                  // z

           colors[3*id]     = 200;                                // r
           colors[3*id+1]   = 200;                                // g
           colors[3*id+2]   = 200;                                // b
         }
         break;
       case 1:
        // urg.get_distance_intensity(dist, intens);
        //TODO:
         cout << "have to been done" << endl;
         break;
       case 2:
        // urg.get_multiecho(dist);
        //TODO:
         cout << "have to been done" << endl;

         break;
       case 3:
         //urg.get_multiecho_intensity(dist, intens);
          // TODO:
         cout << "have to been done" << endl;

       break;
       }


      cloud->setCoords(cloud, cloudsize, 3);
      cloud->setColors(colors, cloudsize, 3);
      viewer3D->update();
     }
  }
};


int main(int argc, char* argv[])
{
  // Sensor initialization
  unsigned int     rays       = hokuyo->getNumberOfRays();
  double  angularRes = hokuyo->getAngularRes();
  double  minPhi     = hokuyo->getStartAngle();

  cout << "rays: " << rays << endl;
  cout << "angularRes: " << angularRes << endl;
  cout << "minPhi: " << minPhi << endl;

/* Init */
  // Load parameters
  init(argc, argv);


  /* Initialization of 3D viewer */
     cloudsize = rays;
     VtkCloud* cloud;

     // size = data_size_Scanner * 3_echos * 3_koordinates
     data   = new double[rays * 3];
     colors = new unsigned char[rays * 3];
     viewer3D->showAxes(true);
     // define update timer
     vtkSmartPointer<vtkTimerCallback> cb =  vtkSmartPointer<vtkTimerCallback>::New();
     vtkSmartPointer<vtkRenderWindowInteractor> interactor = viewer3D->getWindowInteractor();
     interactor->AddObserver(vtkCommand::TimerEvent, cb);
     interactor->CreateRepeatingTimer(30);


     /* Start showing */
           viewer3D->addCloud(cloud);
           viewer3D->startRendering();


/* Shutdown */
   delete hokuyo;
   delete viewer3D;
   delete [] data;
   delete [] colors;
   delete cloud;
   return 0;

}
