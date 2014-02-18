/**
* @file   ${fileShow.h}
* @author Rainer Koch
* @date   ${2014/18/02}
*
* Load different file clouds to filter reflections
*/


void createCube(int cubesize);

/*
 * read from x,y file
 * */
int load_xy_file(char* filename);
/*
 * read from x,y,z file
 * */
int load_xyz_file(char* filename);

/*
 * filter to detect reflections
 * */

void filter(double* distance, unsigned char* colors, double* intensity);


