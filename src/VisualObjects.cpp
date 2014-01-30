/*
 * VisualObjects.cpp
 *
 *  Created on: 30.01.2014
 *      Author: rainer
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <stdlib.h>
#include <math.h>

#include "VisualObjects.h"

#include "obgraphic/Obvious3D.h"
#include <vector>

#define _USE_MATH_DEFINES

using namespace obvious;



VtkCloud* createCupe(int cubesize);
{
  VtkCloud* cloud_cube;
  cloud_cube = new VtkCloud();
//  double* data_cube;
//
//  unsigned char* colors_cube;

  double* data_cube   = new double[3*cubesize*cubesize*cubesize];
  unsigned char* color_cube = new unsigned char[3*cubesize*cubesize*cubesize];


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

        colors_cube[3*n] = 255;
        colors_cube[3*n+1] = 215;
        colors_cube[3*n+2] = 0;
     //   cout << data_cube[3*n] << " " <<  data_cube[3*n+1] << " " << data_cube[3*n+2] << endl;
        n++;
      }
    }
  }

  cloud_cube->setCoords(data_cube, (cubesize*cubesize*cubesize), 3);
  cloud_cube->setColors(color_cube, (cubesize*cubesize*cubesize), 3)

    return cloud_cube;
}


