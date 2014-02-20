/**
* @file   ${fileShow.h}
* @author Rainer Koch
* @date   ${2014/18/02}
*
* Load different file clouds to filter reflections
*/


void createCube(int cubesize);

/*
 * read from txt file
 * x, y
 * */
int load_xy_file(char* filename);
/*
 * read from txt file
 * x, y, intensity
 * */
int load_xyi_file(char* filename);
/*
 * read from txt file
 * 3 echos of x, y
 * */
int load_3xy_file(char* filename);
/*
 * read from txt file
 * 3 echos of x, y, intensity
 * */
int load_xyi_file(char* filename);

/*
 * read from txt file
 * x, y, z
 * */
int load_xyz_file(char* filename);
/*
 * read from rxp file
 * x, y, z, reflectance, amplitude, deviation
 * */
int load_rxp_file(char* filename);

/*
 * filter to detect reflections
 * */
void filter(double* distance, unsigned char* colors, float* intensity);
/*
 * max. limits for data
 */
void boundryBox();

/*
 * reduce shown points
 */
void cbPointIncrease();
void cbPointReduce();
void cbPointReducerEnable();


void initShowCloud();
void createShowCloud();
