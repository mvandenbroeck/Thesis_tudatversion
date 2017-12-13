#include "algorithm"
#include "constants.h"
#include "taylorSeriesIntegratorGenMap.h"
#include "spacecraft.h"
#include "integrationSettings.h"
#include "lagrangeInterpolator.h"
#include "stepSizeControlTSI.h"
#include "problemRecurrenceRelationsGen.h"
#include "discreteForceModelGen.h"

//?
#include <iostream>

#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

#include "boost/make_shared.hpp"

//typedef Eigen::Matrix< double, 8, 1> Vector8d;

TaylorSeriesIntegrator::TaylorSeriesIntegrator(
        bool useCSI,
        bool useCorrectGamma,
        ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer,
        DiscreteForceModelPointer discreteForceModelPointer,
        SpacecraftPointer spacecraftpointer,
        IntegrationSettingsPointer integrationSettingsPointer,
        ConstantsPointer constantsPointer,
        StepSizeControlTSIPointer stepSizeControlTSIPointer):
        useCSI_( useCSI ),
        useCorrectGamma_( useCorrectGamma ),
        problemRecurrenceRelationsPointer_( problemRecurrenceRelationsPointer ),
        discreteForceModelPointer_( discreteForceModelPointer ),
        spacecraftPointer_( spacecraftpointer ),
        integrationSettingsPointer_( integrationSettingsPointer ),
        constantsPointer_( constantsPointer ),
        stepSizeControlTSIPointer_( stepSizeControlTSIPointer){}

TaylorSeriesIntegrator::~TaylorSeriesIntegrator()
{

}


void TaylorSeriesIntegrator::integrate()
{

    int numberOfStateVariables = spacecraftPointer_->getInitialState().rows();
    Eigen::MatrixXd initialState( numberOfStateVariables, 1 );
    initialState = spacecraftPointer_->getInitialState();

    // Initiate values for the while loop
    double currentTime = integrationSettingsPointer_->getInitialTime();
    Eigen::MatrixXd currentState(  numberOfStateVariables, 1 );
    currentState = initialState;
    double currentStepSize = integrationSettingsPointer_->getInitialStepSize();
    // Save currentStepSize in integrationSettingsPointer
    integrationSettingsPointer_->setStepSize( currentStepSize );

    // Declare stateMap and thrustAccelerationMap to store output
    std::map< double, Eigen::MatrixXd > computedThrustAccelerationMap;
    std::map< double, Eigen::MatrixXd > stateMap;
    std::map< double, Eigen::MatrixXd > fullStateMap;
    stateMap[ currentTime ] = initialState;
    fullStateMap[ currentTime ] = initialState;

    //??
//    std::cout << "stateMap[ " << currentTime << " ] =\n" << stateMap[ currentTime ] << std::endl;

    // Define thrustForceMatrix for in loop use
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix;
    thrustForceMatrix = spacecraftPointer_->getThrustForceMatrix();

    // Define variables for in loop use
    double finalTime = integrationSettingsPointer_->getFinalTime();
    int order = integrationSettingsPointer_->getOrderOfTaylorSeries();
    int numberOfVariables = spacecraftPointer_->getInitialReducedState().size();
    double minimumStepSize = integrationSettingsPointer_->getMinimumStepSize();
    double maximumStepSize = integrationSettingsPointer_->getMaximumStepSize();
    int refineFactor = integrationSettingsPointer_->getRefineFactor();

    if ( useCSI_ == true )
    {
        // Define time and thrust values in std::vector
        std::vector< double > timePointsToBeInterpolated( thrustForceMatrix.rows() );
        std::vector< double > thrustPointsToBeInterpolatedX( thrustForceMatrix.rows() );
        std::vector< double > thrustPointsToBeInterpolatedY( thrustForceMatrix.rows() );
        std::vector< double > thrustPointsToBeInterpolatedZ( thrustForceMatrix.rows() );
        std::vector< double > resultantThrustPointsToBeInterpolated( thrustForceMatrix.rows() );

        for ( int i = 0; i < thrustForceMatrix.rows(); i++ )
        {
            timePointsToBeInterpolated[ i ] = thrustForceMatrix( i, 0 );
            thrustPointsToBeInterpolatedX[ i ] = thrustForceMatrix( i, 1 );
            thrustPointsToBeInterpolatedY[ i ] = thrustForceMatrix( i, 2 );
            thrustPointsToBeInterpolatedZ[ i ] = thrustForceMatrix( i, 3 );
            resultantThrustPointsToBeInterpolated[ i ] = ( ( thrustForceMatrix.rightCols( 3 ) ).row( i ) ).norm();
        }

        // Calculate the four coefficients of all splines (number of splines = 1 - number of nodes)
        // of the cubic spline interpolation.
        tudat::interpolators::CubicSplineInterpolatorDoublePointer thrustInVelocicityFrameXcsiPointer =
                boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                    timePointsToBeInterpolated,
                    thrustPointsToBeInterpolatedX );
        tudat::interpolators::CubicSplineInterpolatorDoublePointer thrustInVelocicityFrameYcsiPointer =
                boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                    timePointsToBeInterpolated,
                    thrustPointsToBeInterpolatedY );
        tudat::interpolators::CubicSplineInterpolatorDoublePointer thrustInVelocicityFrameZcsiPointer =
                boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                    timePointsToBeInterpolated,
                    thrustPointsToBeInterpolatedZ );
        tudat::interpolators::CubicSplineInterpolatorDoublePointer resultantThrustForceCSIpointer =
                boost::make_shared< tudat::interpolators::CubicSplineInterpolatorDouble >(
                    timePointsToBeInterpolated,
                    resultantThrustPointsToBeInterpolated );

        std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX,
                                             CSIcoefficientsOfThrustForceInVelocityFrameY,
                                             CSIcoefficientsOfThrustForceInVelocityFrameZ,
                                             CSIcoefficientsOfResultantThrustForce;
        CSIcoefficientsOfThrustForceInVelocityFrameX = thrustInVelocicityFrameXcsiPointer->GetCoefficients();
        CSIcoefficientsOfThrustForceInVelocityFrameY = thrustInVelocicityFrameYcsiPointer->GetCoefficients();
        CSIcoefficientsOfThrustForceInVelocityFrameZ = thrustInVelocicityFrameZcsiPointer->GetCoefficients();
        CSIcoefficientsOfResultantThrustForce = resultantThrustForceCSIpointer->GetCoefficients();


        // Store coefficients in spacecraft vector
        spacecraftPointer_->setCSIcoefficientsOfThrustForceInVelocityFrameX( CSIcoefficientsOfThrustForceInVelocityFrameX );
        spacecraftPointer_->setCSIcoefficientsOfThrustForceInVelocityFrameY( CSIcoefficientsOfThrustForceInVelocityFrameY );
        spacecraftPointer_->setCSIcoefficientsOfThrustForceInVelocityFrameZ( CSIcoefficientsOfThrustForceInVelocityFrameZ );
        spacecraftPointer_->setCSIcoefficientsOfResultantThrustForce( CSIcoefficientsOfResultantThrustForce );
    }

    Eigen::MatrixXd dummy_variable;

    // Initiate coefficientMatrix and store in spacecraft pointer
    Eigen::MatrixXd coefficientMatrix( order + 1, numberOfStateVariables );
    if ( useCorrectGamma_ == true )
    {
        coefficientMatrix.fill( 0.0 ); // Only for the first step all the derivatives are zero.
        spacecraftPointer_->setCurrentCoefficientMatrix( coefficientMatrix );
        // the first row can be filled with the state vector, but this is not needed as this row is not used

    }

    // Declarations for variables used in time loop
    Eigen::MatrixXd currentThrustAccelerationVector( 1, 3 );
    double previousTime;
    Eigen::MatrixXd  nextStateVector(  numberOfStateVariables, 1 );
    Eigen::MatrixXd interTime( order + 1, 1 );
    Eigen::MatrixXd coefficientMatrixInLoop( order + 1, 1 );


    // Advance the state
    while( currentTime < finalTime ) // while the current time is smaller than the final time
    {
        //?? Display currentTime
//        std::cout << "currentTime = " << currentTime << std::endl;


        /// Use discrete force model to update the current force,
        /// force derivatives, acceleration and acceleration derivatives in the spacecraft body
        // Update variables in spacecraft body
        discreteForceModelPointer_->updateCurrentForcesAndAccelerationsForTSI( currentTime ); // Writes output in spacecraftPointer_
        // is not done for the finalTime. So at finalTime the current forces and accelerations are not calculated

        // Get current thrust accelerationvector from spacecraft pointer. Must be computed before next step size is calculated
        currentThrustAccelerationVector = spacecraftPointer_->getCurrentAcceleration();

        // Add current thrust acceleration to the map
        computedThrustAccelerationMap[ currentTime ] = currentThrustAccelerationVector;


        /// The general algorithm starts here (problem independent)

        // Store currentTime for the refinement
        previousTime = currentTime;

        // Advance the time
        currentTime = currentTime + currentStepSize;

        if( currentTime > finalTime ) // if the next time is bigger than the final time
        {
            currentStepSize = currentStepSize - ( currentTime - finalTime ); // get the correct step size
            currentTime = finalTime; // get the correct time
        }

        // Calculate next state and next time and matrix of derivatives. The output will be saved in
        // the spacecaftPointer_
        problemRecurrenceRelationsPointer_->computeNextState( currentStepSize );

        //! Get output of problemRecurrenceRelations
        // Get current state vector from spacecraft pointer
        nextStateVector = spacecraftPointer_->getCurrentState();
        // Get current coefficients of derivative matrix from spacecraft pointer.
        coefficientMatrix = spacecraftPointer_->getCurrentCoefficientMatrix();

        // Save the next state in the map along with the updated currentTime
        //??
//        if ( currentTime == 0.0 )
//        {
//            std::cout << "Label 1" << std::endl;
//            std::cout << "stateMap[ " << currentTime << " ] =\n" << stateMap[ currentTime ] << std::endl;
//        }
        stateMap[ currentTime ] = nextStateVector;

        // Define new matrix with the coefficients of the Taylor series
        Eigen::MatrixXd coefficientMatrixReduced( coefficientMatrix.rows(), numberOfVariables );
        coefficientMatrixReduced = coefficientMatrix.leftCols( numberOfVariables );

        /// Next time step
        // Define next step size
        double temp_currentStepSize = currentStepSize;
        currentStepSize = stepSizeControlTSIPointer_->determineNextStepSize( coefficientMatrixReduced, temp_currentStepSize );


        // Correct step size in case step size is outside limits
        if ( currentStepSize < minimumStepSize )
        {
            currentStepSize = minimumStepSize;
        }
        else if ( currentStepSize > maximumStepSize )
        {
            currentStepSize = maximumStepSize;
        }

        // Save next step size in integrationSettingsPointer
        integrationSettingsPointer_->setStepSize( currentStepSize );

        // Refinement
        if ( refineFactor > 0 ) // refineFactor is defined as an int and defines amount of points in between
            // previousTime and currentTime
        {
            double referenceStepSize = ( currentTime - previousTime ) / ( refineFactor + 1 );

            // Compute the (refineFactor) additional steps in between previous and current time
            Eigen::VectorXd referenceTime( refineFactor );
            referenceTime.setLinSpaced( refineFactor, previousTime + referenceStepSize, currentTime - referenceStepSize );

            Eigen::MatrixXd interStatesMatrix( refineFactor, numberOfVariables );
            //interStatesMatrix.fill( 0.0 );

            for ( int n = 0; n < refineFactor; n++ )
            {
                for ( int orderInLoop = 0; orderInLoop <= order; orderInLoop++ )
                {
                    interTime( orderInLoop, 0 ) = std::pow( referenceTime( n ) - previousTime , orderInLoop );
                }

                for ( int n2 = 0; n2 < numberOfVariables; n2++ )
                {
                    coefficientMatrixInLoop = coefficientMatrix.col( n2 );
                    dummy_variable = coefficientMatrixInLoop.adjoint() * interTime;
                    interStatesMatrix( n, n2 ) = dummy_variable( 0 );
                }
            }

            Eigen::MatrixXd newInterStateMatrix( refineFactor, numberOfStateVariables );
            newInterStateMatrix << interStatesMatrix, Eigen::MatrixXd::Zero( refineFactor, ( numberOfStateVariables - numberOfVariables ) );

//            // Add refined states to state matrix
//            Eigen::MatrixXd temp_stateMatrix( stateMatrix.rows() + refineFactor, stateMatrix.cols() );
//            temp_stateMatrix << stateMatrix.topRows( stateMatrix.rows() - 1 ),
//                                newInterStateMatrix,
//                                stateMatrix.bottomRows( 1 );

//            stateMatrix = temp_stateMatrix;

//            // Add refined times to state matrix
//            Eigen::MatrixXd temp_Time( Time.rows() + refineFactor, Time.cols() );
//            temp_Time << Time.topRows( Time.rows() - 1 ),
//                         referenceTime,
//                         Time.bottomRows( 1 );

//            Time = temp_Time;

            // First add refined steps, then add the calculated step.
            for ( int i = 0; i < refineFactor; i++ )
            {
                fullStateMap[ referenceTime( i ) ] << newInterStateMatrix.row( i );
            }
            // Finally add the currentTime and currentState
            fullStateMap[ currentTime ] = currentState;
        }

    }

//    //??
//    std::cout << "stateMap[ 0 ] =\n" << stateMap[ 0.0 ] << std::endl;

    // Store output
    setStateMap( stateMap );
    setComputedThrustAccelerationMap( computedThrustAccelerationMap );
    setFullStateMap( fullStateMap ); // inclusive the refined states
}
