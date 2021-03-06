/*
 * LaserDevice.h
 *
 *  Created on: 11.03.2013
 *      Author: Rainer
 */

#ifndef LASER2D_H_
#define LASER2D_H_

/**
 * @namespace obvious
 */
namespace laserDev
{
/**
 * @class LaserDevice
 */
class LaserDevice
{
public:
  /**
   * Standard Constructor
   */
  LaserDevice();

  /**
   * Default Destructor
   */
  virtual   ~LaserDevice(void)
  {
    delete[] _ranges;
    delete[] _intensities;
    delete[] _angles;
    delete[] _mask;
  };

  /**
   * Function to return start Angle
   * @return  angle in rad
   */
  double UTM30LXEW::getStartAngle();

  /**
   * Function to return stop Angle
   * @return  angle in rad
   */
  double UTM30LXEW::getStopAngle();

  unsigned int getNumberOfRays();

  double getAngularRes(void) {return _nrOfRays);

  double* getRanges();

  void schedule();

  /**
   * Function to grab new data
   * @return  TRUE if success
   */
  virtual bool grab(void)           { return(true); }

  /**
   * Function to return Distances
   * @return  Distance in meters
   */
//  double*   getDistance(void) const    { return _ranges; }

  /**
   * Function to return coords in 2d (x1,y1,x2,y2, ...)
   * @return  coords2D
   */
  double*   calculateCoords2D(void) const    { return _coords2D; }

  /**
   * Function to return coords in 3d (x1,y1,z1, ...)
   * @return  coords3D
   *
   * Please note that the paramter z1 is set to zero.
   */
  double*   estimateCoords2D(void) const     { return _coords3D; }

  /**
   * Function to return normals
   * @return  normals in 2d vector
   */
//  double*   getNormals(void) const        { return _normals; }

  /**
   * Function to return intensities
   * @return  Intensity value
   */
  double*   calculateIntensities(void) const   { return _intensities; }

  /**
   * Function to return angles
   * @return  angles in rad
   */
  double*   calculateAngles(void) const        { return _angles; }

  /**
   * Function to return mask of valid points
   * @return  mask with valid points (TRUE)
   */
//  bool*     getMask(void)                { return _mask; }

protected:
  /**
   * Function to estimate ranges in scan
   */
  virtual void calculateRanges(void) { }

  /**
   * Function to estimate intensities in scan
   */
  virtual void calculateIntensities(void) { }

  /**
   * Function to estimate single angles for every ray
   */
  virtual void calculateAngles(void) { }

  /**
   * Function to estimate 2D coords
   */
  virtual void estimateCoords2D(void) { }

  /**
   * Function to estimate 2D coords
   */
  virtual void estimateCoords3D(void) { }

  /**
   * Function to estimate normals
   */
  //virtual void estimateNormals(void) { }

  /**
   * Function to estimate mask
   */
//  virtual void estimateMask() { }

  /*~~~~~~~~~~~~~~~~~~~~~~~~~~ Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  double*   _ranges;            //!< Distance in meters
  double*   _intensities;       //!< Intensities
  double*   _coords2D;          //!< 2D coords
  //double*   _coords3D;          //!< 3D coords
  //double*   _normals;           //!< normals
  double*   _angles;            //!< Angles in rad
  bool*     _mask;              //!< mask for valid or invalid points
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~ Configuration ~~~~~~~~~~~~~~~~~~~~~~~~~*/

 double      _nrOfRays;
  double            _angRes;

};

};  //namespace


#endif /* LASERDEVICE_H_ */
