#include "frameTransformation.h"
#include "Eigen/Core"
#include "Eigen/Geometry"

//??
#include <limits>
#include <iomanip>
#include <iostream>

FrameTransformation::FrameTransformation() {}

Eigen::Matrix< double, 1, 3 > FrameTransformation::velocityFrameToUSMFrame(
        Eigen::Matrix< double, 1, 3> thrustInVelocityFrame,
        double flightPathAngle )
{
    Eigen::Matrix< double, 1, 3 > thrustInUSMFrame;
    thrustInUSMFrame( 0, 0 ) = - thrustInVelocityFrame( 0, 1 ) * std::cos( flightPathAngle ) +
            thrustInVelocityFrame( 0, 0 ) * std::sin( flightPathAngle );
    thrustInUSMFrame( 0, 1 ) = thrustInVelocityFrame( 0, 1 ) * std::sin( flightPathAngle ) +
            thrustInVelocityFrame( 0, 0 ) * std::cos( flightPathAngle );
    thrustInUSMFrame( 0, 2 ) = thrustInVelocityFrame( 0, 2 );

    return thrustInUSMFrame;
}

Eigen::Matrix< double, 1, 3 > FrameTransformation::velocityFrameToCartFrame(Eigen::Matrix< double, 1, 3 > thrustInVelocityFrame,
        Eigen::Matrix<double, 6, 1 > vehicleStateCartesian )
{
    Eigen::Vector3d vehicleVelocity, vehicleRadius;
    vehicleRadius = vehicleStateCartesian.head( 3 );
    vehicleVelocity = vehicleStateCartesian.tail( 3 );

    Eigen::Vector3d unitT = vehicleVelocity / vehicleVelocity.norm( );
    if ( vehicleRadius.cross( vehicleVelocity ).norm( ) == 0.0 )
    {
        std::string errorMessage = "Division by zero: radius and velocity are in the same direction.";
        throw std::runtime_error( errorMessage );
    }

    Eigen::Vector3d unitW = ( vehicleRadius.cross( vehicleVelocity ) ).normalized( );

    Eigen::Vector3d unitN = ( unitW.cross( unitT ) ).normalized( );

    Eigen::Matrix3d transformationMatrix;
    transformationMatrix << unitT( 0 ), unitN( 0 ), unitW( 0 ),
                            unitT( 1 ), unitN( 1 ), unitW( 1 ),
                            unitT( 2 ), unitN( 2 ), unitW( 2 );

    //??
//    std::cout << "transformationMatrix =\n" << std::setprecision( std::numeric_limits< double >::max_digits10 )
//              << transformationMatrix << std::endl;

    //??
//    std::cout << "Thrust in VF =\n" << std::setprecision( std::numeric_limits< double >::max_digits10 )
//              << thrustInVelocityFrame << std::endl;

    Eigen::Matrix< double, 1, 3 > thrustInCartFrame =
             ( transformationMatrix * thrustInVelocityFrame.transpose() ).transpose();

    //??
//    std::cout << "Thrust in Cart frame =\n" << std::setprecision( std::numeric_limits< double >::max_digits10 )
//              << thrustInCartFrame << std::endl;

    return thrustInCartFrame;
}

FrameTransformation::~FrameTransformation()
{

}

