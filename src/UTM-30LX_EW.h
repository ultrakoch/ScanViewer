/**
* @file   ${UTM-30LX_EW.h}
* @author Rainer Koch
* @date   ${2014/01/09}
*
* Class for Hokuyo UTM-30LX_EW Laser scanner
*/



#ifndef UTM30LXEW_H_
#define UTM30LXEW_H_

#include "Lidar.h"
#include <Urg_driver.h>



/**
 * @namespace obvious
 */
namespace laserDev
{
  /**
   * @class LaserDevice
   */
  class UTM30LXEW
  {
  public:
    /**
     * Standard Constructor
     */
    UTM30LXEW();

    /**
     * Default Destructor
     */
    virtual   ~UTM30LXEW();

    double getStartAngle();

    double getStopAngle();

    unsigned int getNumberOfRays();

    double getAngularRes(void);

    double* getRanges();

    double* getCoords2D();

    /**
     * Function to grab new data
     * @return  TRUE if success
     */
    virtual bool      grab(void);

    void schedule();

  private:

    /**
     * Function to estimate ranges in scan
     */
    void calculateRanges(void);

    /**
     * Function to estimate intensities in scan
     */
    void calculateIntensities(void);

    /**
     * Function to estimate single angles for every ray
     */
    void calculateAngles(void);

    /**
     * Function to estimate 2D coords
     */
    void calculateCoords2D(void);



    qrk::Urg_driver  _laser;

    unsigned int _nrOfRays;

    double*   _ranges;            //!< Distance in meters
    double*   _intensities;       //!< Intensities
    double*   _coords2D;          //!< 2D coords
    double*   _normals;           //!< normals
    double*   _angles;            //!< Angles in rad
    bool*     _mask;              //!< mask for valid or invalid points
    std::vector<long> dist;
    std::vector<long> _distBuffer;

    std::vector<unsigned short> intens;


  };

};  //namespace


#endif /* UTM30LXEW_H_ */
