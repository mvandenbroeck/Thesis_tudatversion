#ifndef SPACECRAFT_H
#define SPACECRAFT_H

#include <array>
#include <vector>
#include <Eigen/Core>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include "constants.h"

#include <boost/shared_ptr.hpp>



class Spacecraft
{
public:
    Spacecraft(
            double initialMass,
            Vector7d initialReducedState,
            double specificImpulse,
            ConstantsPointer constantsPointer);
    double getMass() {return mass_;}
    double getSpecificImpulse() {return specificImpulse_;}
    // Add a getTotalState function that gives the state at all timesteps
    Vector11d getInitialState() {return initialState_;}
    Vector7d getInitialReducedState() {return initialReducedState_;}
    Eigen::MatrixXd getCurrentAcceleration() {return currentAcceleration_;}
    Eigen::MatrixXd getCurrentAccelerationDerivative() {return currentAccelerationDerivative_;}
    Eigen::MatrixXd getCurrentResultantThrustForceDerivatives() {return currentResultantThrustForceDerivatives_;}
    Eigen::MatrixXd getCurrentState() {return currentState_;}
    Eigen::MatrixXd getCurrentCoefficientMatrix() {return currentCoefficientMatrix_;}
    Eigen::Matrix< double, 1, 3 > getCurrentThrustVector() {return currentThrustVector_;}
    double getCurrentTime() {return currentTime_;}
    std::vector< std::vector< double > > getCSIcoefficientsOfThrustForceInVelocityFrameX()
                                {return CSIcoefficientsOfThrustForceInVelocityFrameX_;}
    std::vector< std::vector< double > > getCSIcoefficientsOfThrustForceInVelocityFrameY()
                                {return CSIcoefficientsOfThrustForceInVelocityFrameY_;}
    std::vector< std::vector< double > > getCSIcoefficientsOfThrustForceInVelocityFrameZ()
                                {return CSIcoefficientsOfThrustForceInVelocityFrameZ_;}
    std::vector< std::vector< double > > getCSIcoefficientsOfResultantThrustForce()
                                {return CSIcoefficientsOfResultantThrustForce_;}
    double getCurrentResultantThrustForce()
                                {return currentResultantThrustForce_;}

    void setMass(double mass) {mass_ = mass;}
    void setSpecificImpulse(double specificImpulse) {specificImpulse_ = specificImpulse;}
    void loadThrustForce(std::string fileName);
    Eigen::Matrix< double, Eigen::Dynamic, 4 > getThrustForceMatrix() {return thrustForceMatrix_;}

    void setCurrentAcceleration( Eigen::MatrixXd currentAcceleration)
                                {currentAcceleration_ = currentAcceleration;}
    void setCurrentAccelerationDerivative( Eigen::MatrixXd currentAccelerationDerivative )
                                {currentAccelerationDerivative_ = currentAccelerationDerivative;}
    void setCurrentResultantThrustForceDerivatives( Eigen::MatrixXd currentResultantThrustForceDerivatives ) // dim max( order, 3) x 1
                                {currentResultantThrustForceDerivatives_ = currentResultantThrustForceDerivatives;}
    void setCurrentState( Vector11d currentState )
                                {currentState_ = currentState;}
    void setCurrentCoefficientMatrix( Eigen::MatrixXd currentCoefficientMatrix )
                                {currentCoefficientMatrix_ = currentCoefficientMatrix;}
    void setCurrentFlightPathAngle( double currentFlightPathAngle )
                                {currentFlightPathAngle_ = currentFlightPathAngle;}
    void setCurrentTime( double currentTime ) {currentTime_ = currentTime;}
    void setCSIcoefficientsOfThrustForceInVelocityFrameX( std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX )
                                {CSIcoefficientsOfThrustForceInVelocityFrameX_ = CSIcoefficientsOfThrustForceInVelocityFrameX;}
    void setCSIcoefficientsOfThrustForceInVelocityFrameY( std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameY )
                                {CSIcoefficientsOfThrustForceInVelocityFrameY_ = CSIcoefficientsOfThrustForceInVelocityFrameY;}
    void setCSIcoefficientsOfThrustForceInVelocityFrameZ( std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameZ )
                                {CSIcoefficientsOfThrustForceInVelocityFrameZ_ = CSIcoefficientsOfThrustForceInVelocityFrameZ;}
    void setCSIcoefficientsOfResultantThrustForce( std::vector< std::vector< double > > CSIcoefficientsOfResultantThrustForce )
                                {CSIcoefficientsOfResultantThrustForce_ = CSIcoefficientsOfResultantThrustForce;}
    void setCurrentResultantThrustForce( double currentResultantThrustForce )
                                {currentResultantThrustForce_ = currentResultantThrustForce;}

private:
    double mass_;
    double specificImpulse_;
    ConstantsPointer constantsPointer_;
    Vector7d initialReducedState_;

    Vector11d initialState_;
    Vector11d currentState_;

    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix_;

    Eigen::MatrixXd currentAcceleration_;
    Eigen::MatrixXd currentAccelerationDerivative_;
    Eigen::MatrixXd currentResultantThrustForceDerivatives_;
    Eigen::MatrixXd currentCoefficientMatrix_;
    Eigen::Matrix< double, 1, 3 > currentThrustVector_;
    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX_;
    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameY_;
    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameZ_;
    std::vector< std::vector< double > > CSIcoefficientsOfResultantThrustForce_;
    double currentResultantThrustForce_;
    double currentFlightPathAngle_;
    double currentTime_;


};

typedef boost::shared_ptr< Spacecraft > SpacecraftPointer;

#endif // SPACECRAFT_H
