/**
* @file   ${hokuyoShow.h}
* @author Rainer Koch
* @date   ${2013/12/06}
*
* Stream Hokuyo to file
*/

#include "Urg_driver.h"

void save_raw_dist(std::vector<long>& distance);
void save_xy_dist(std::vector<long>& distance);
void save_raw_multidist(std::vector<long>& distance);
void save_xy_multidist(std::vector<long>& distance);
void save_raw_intens(std::vector<long>& distance, std::vector<unsigned short>& intensity);
void save_xy_intens(std::vector<long>& distance, std::vector<unsigned short>& intensity);
void save_raw_mulitintens(std::vector<long>& distance, std::vector<unsigned short>& intensity);
void save_xy_mulitintens(std::vector<long>& distance, std::vector<unsigned short>& intensity);

void createCupe(int cubesize);

/* read from x,y file
 *
 * */
int load_xy_file(char* filename);

void startLaser();
void stopLaser();

void cbSave();
void cbMeasuretype();



