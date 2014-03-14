/**
* @file   ${fileShow.h}
* @author Rainer Koch
* @date   ${2014/18/02}
*
* Load different file clouds to filter reflections
*/

void testfkt();
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
void filter(double* distance, unsigned char* colors, double* intensity, int cloudsize);
void deltectMirroredPoints(double* cloud, int size, double origin[3], double axis1[3], double axis2[3]);
double* findPlane(double* cloud, int size, double* planeCenterVektor);

/*
 * Functions used in filter to mark the points depending on the location
 * */
void setAsMirror(int n, double* testpoint);
void setAsObject(int n, double* testpoint);
void setAsMirror(int n, double* testpoint);

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
/*
 * mark the data point with a color depending on his intensity
 */
void cbIntensityColor();

void initShowCloud();
void createShowCloud();


