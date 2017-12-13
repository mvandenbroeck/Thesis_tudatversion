#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Core>

#include <tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h>



#include "USMstateDerivativeModel.h"
#include "constants.h"
#include "spacecraftGen.h"
#include "frameTransformation.h"


USMstateDerivativeModel::USMstateDerivativeModel(
        SpacecraftGenPointer spacecraftGenPointer,
        ConstantsPointer constantsPointer ) :
    spacecraftGenPointer_( spacecraftGenPointer ),
    constantsPointer_( constantsPointer ) {}

void USMstateDerivativeModel::computeCurrentAccelerations( const double &currentTime )
{
    /// Extract variables from classes
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix;
    thrustForceMatrix = spacecraftGenPointer_->getThrustForceMatrix();

    /// Specific relations for the specific low-thrust representation method
    // Find index in thrust profile that corresponds with the current time
    int index = tudat::basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                thrustForceMatrix.col( 0 ), currentTime);
    double Time_0 = thrustForceMatrix( index, 0 ); // calculate time of left node

    // Get CSI coefficients of resultant and x, y and z components in velocity frame from spacecraft pointer
    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX,
                                         CSIcoefficientsOfThrustForceInVelocityFrameY,
                                         CSIcoefficientsOfThrustForceInVelocityFrameZ,
                                         CSIcoefficientsOfResultantThrustForce;

    CSIcoefficientsOfThrustForceInVelocityFrameX =
            spacecraftGenPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameX(  );
    CSIcoefficientsOfThrustForceInVelocityFrameY =
            spacecraftGenPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameY(  );
    CSIcoefficientsOfThrustForceInVelocityFrameZ =
            spacecraftGenPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameZ(  );
    CSIcoefficientsOfResultantThrustForce =
            spacecraftGenPointer_->getCSIcoefficientsOfResultantThrustForce(  );

    /// Calculate values of resultant and x, y and z components in velocity frame
    double time_temp = currentTime - Time_0;
    double currentThrustForceInVelocityFrameX =
              CSIcoefficientsOfThrustForceInVelocityFrameX[ 0 ][ index ]
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 1 ][ index ] * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ] * time_temp * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameY =
            CSIcoefficientsOfThrustForceInVelocityFrameY[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameZ =
            CSIcoefficientsOfThrustForceInVelocityFrameZ[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentResultantThrustForce = std::sqrt(
            currentThrustForceInVelocityFrameX * currentThrustForceInVelocityFrameX +
            currentThrustForceInVelocityFrameY * currentThrustForceInVelocityFrameY +
            currentThrustForceInVelocityFrameZ * currentThrustForceInVelocityFrameZ );

    /// Calculation of the current thrust acceleration in USM frame components

    // Extract mass from current state
    Eigen::MatrixXd currentState = spacecraftGenPointer_->getCurrentState();
    double currentMass = currentState( 7 );

    // Calculate flight path angle gamma required for transformation from VF to USM frame
//    double ve1 = currentState( 8 );
//    double ve2 = currentState( 9 );
    double C, Rf1, Rf2, e3, eta;
    C = currentState( 0 );
    Rf1 = currentState( 1 );
    Rf2 = currentState( 2 );
    e3 = currentState( 5 );
    eta = currentState( 6 );
    double cosineLambda, sineLambda;
    cosineLambda = ( eta * eta - e3 * e3 ) / ( e3 * e3 + eta * eta );
    sineLambda = ( 2 * e3 * eta ) / ( e3 * e3 + eta * eta );
    double ve1, ve2;
    ve1 = cosineLambda * Rf1 + sineLambda * Rf2;
    ve2 = C - sineLambda * Rf1 + cosineLambda * Rf2;
    double sin_gamma, cos_gamma, gamma;
    sin_gamma = ve1 / std::sqrt( ve1 * ve1 + ve2 * ve2 );
    cos_gamma = ve2 / std::sqrt( ve1 * ve1 + ve2 * ve2 );
    gamma = std::atan2( sin_gamma, cos_gamma );

    // Create frameTransformation object
    FrameTransformation frameTransformation;

    // Create vector of currentThrustForceInVelocityFrame
    Eigen::Matrix< double, 1, 3 > currentThrustForceInVelocityFrame;
    currentThrustForceInVelocityFrame << currentThrustForceInVelocityFrameX,
                                         currentThrustForceInVelocityFrameY,
                                         currentThrustForceInVelocityFrameZ;

    // Compute current thrust force in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustForceInUSMframe =
            frameTransformation.velocityFrameToUSMFrame( currentThrustForceInVelocityFrame, gamma );

    // Compute current thrust acceleration in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustAccelerationInUSMframe =
            currentThrustForceInUSMframe / currentMass;

    // Store current thrust acceleration in USM frame in spacecraft pointer
    spacecraftGenPointer_->setCurrentAcceleration( currentThrustAccelerationInUSMframe );
    spacecraftGenPointer_->setCurrentResultantThrustForce( currentResultantThrustForce  );
}


Eigen::VectorXd USMstateDerivativeModel::computeUnifiedStateModelStateDerivative( const double currentTime, const Eigen::VectorXd& currentState )
{
    double C, Rf1, Rf2, e1, e2, e3, eta;
    C = currentState( 0 );
    Rf1 = currentState( 1 );
    Rf2 = currentState( 2 );
    e1 = currentState( 3 );
    e2 = currentState( 4 );
    e3 = currentState( 5 );
    eta = currentState( 6 );
    // mass is not needed in calculation
//    ve1 = currentState( 8 );
//    ve2 = currentState( 9 );
//    omega3 = currentState( 10 );

    // Compute current accelerations using Discrete Thrust Model
    computeCurrentAccelerations( currentTime );
    Eigen::MatrixXd currentAcceleration = spacecraftGenPointer_->getCurrentAcceleration();
    double ae1, ae2, ae3;
    ae1 = currentAcceleration( 0 );
    ae2 = currentAcceleration( 1 );
    ae3 = currentAcceleration( 2 );

    // Calculate auxiliary variables
    double cosineLambda, sineLambda;
    cosineLambda = ( eta * eta - e3 * e3 ) / ( e3 * e3 + eta * eta );
    sineLambda = ( 2.0 * e3 * eta ) / ( e3 * e3 + eta * eta );
    double ve2;
    ve2 = C - sineLambda * Rf1 + cosineLambda * Rf2;
    double p, k, omega1, omega3;
    p = C / ve2;
    k = ( e1 * e3 - e2 * eta ) / ( e3 * e3 + eta * eta );
    omega1 = ae3 / ve2;
    omega3 = C * ve2 * ve2 / constantsPointer_->centralBodyGravitationalParameter_;

    // Calculate derivative of USM elemenets
    double C_, Rf1_, Rf2_, e1_, e2_, e3_, eta_, m_;
    C_ = -p * ae2;
    Rf1_ = ae1 * cosineLambda - ae2 * ( 1 + p ) * sineLambda - ae3 * k * Rf2 / ve2;
    Rf2_ = ae1 * sineLambda + ae2 * ( 1 + p ) * cosineLambda + ae3 * k * Rf1 / ve2;
    e1_ = 0.5 * ( omega3 * e2 + omega1 * eta );
    e2_ = 0.5 * ( -omega3 * e1 + omega1 * e3 );
    e3_ = 0.5 * ( -omega1 * e2 + omega3 * eta );
    eta_ = 0.5 * ( -omega1 * e1 - omega3 * e3 );
    m_ =  - spacecraftGenPointer_->getCurrentResultantThrustForce() /
            ( constantsPointer_->standardGravity_ * spacecraftGenPointer_->getSpecificImpulse() );
//    ve1_ = omega3 * ve2 - omega3 * C;
//    ve2_ = -omega3 * ve1;
//    omega3_ = 2.0 * omega3 * ve2_ / ve2;

    // Save derivatives in outputvector
    Eigen::VectorXd Derivative( 8 );
    Derivative << C_, Rf1_, Rf2_, e1_, e2_, e3_, eta_, m_;
    return Derivative;
}
