#ifndef CONSTANTS
#define CONSTANTS

#include "Eigen/Core"
#include <boost/shared_ptr.hpp>

// It is better to use a constants class instead of globals because if the constants have to be changed, then the whole
// project has to be built again, which takes a lot of time.

// using thrustForceMatrix = std::vector<std::array<double,3>>;
typedef Eigen::Matrix< double, 11, 1> Vector11d;
typedef Eigen::Matrix< double, 8, 1> Vector8d;
typedef Eigen::Matrix< double, 7, 1> Vector7d;
//typedef Eigen::Matrix< double, Eigen::Dynamic, 3> MatrixX3d;




class Constants
{
public:
    Constants( const double centralBodyGravitationalParameter,
               const double standardGravity ):
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        standardGravity_( standardGravity ){}

    const double centralBodyGravitationalParameter_ = 1.32712440018e20;  // [m^3/s^2] // 1.32712440018e20
    const double standardGravity_ = 9.80665; // [m/s^2] // 9.80665
    const double astronomicalUnit_ = 1.49597870691e11; // [m]
};

typedef boost::shared_ptr< Constants > ConstantsPointer;

#endif // CONSTANTS

