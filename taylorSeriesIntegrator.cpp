#include "algorithm"
#include "constants.h"
#include "taylorSeriesIntegrator.h"
#include "spacecraft.h"
#include "integrationSettings.h"
#include "lagrangeInterpolator.h"
#include "stepSizeControlTSI.h"
#include "problemRecurrenceRelations.h"
#include "discreteForceModel.h"

#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

#include "boost/make_shared.hpp"

//typedef Eigen::Matrix< double, 8, 1> Vector8d;

TaylorSeriesIntegrator::TaylorSeriesIntegrator(
        ProblemRecurrenceRelationsPointer problemRecurrenceRelationsPointer,
        SpacecraftPointer spacecraftpointer,
        IntegrationSettingsPointer integrationSettingsPointer,
        ConstantsPointer constantsPointer,
        StepSizeControlTSIPointer stepSizeControlTSIPointer):
        problemRecurrenceRelationsPointer_( problemRecurrenceRelationsPointer ),
        spacecraftPointer_( spacecraftpointer ),
        integrationSettingsPointer_( integrationSettingsPointer ),
        constantsPointer_( constantsPointer ),
        stepSizeControlTSIPointer_( stepSizeControlTSIPointer){}

TaylorSeriesIntegrator::~TaylorSeriesIntegrator()
{

}

void TaylorSeriesIntegrator::integrate()
{
    // Setting the first time and state outputs
    Eigen::MatrixXd Time( 1, 1 );
    Time( 0, 0 ) = integrationSettingsPointer_->getInitialTime();

    int numberOfStateVariables = spacecraftPointer_->getInitialState().rows();
    Eigen::MatrixXd initialState( numberOfStateVariables, 1 );
    initialState = spacecraftPointer_->getInitialState();

    Eigen::MatrixXd stateMatrix( 1, numberOfStateVariables );
    stateMatrix.row( 0 ) = initialState.adjoint();

    // Declare computedThrustAccelerationMatrix
    Eigen::MatrixXd computedThrustAccelerationMatrix( 1, 3 );

    // Initiate values for the while loop
    double currentTime = integrationSettingsPointer_->getInitialTime();
    Eigen::MatrixXd currentState(  numberOfStateVariables, 1 );
    currentState = initialState;
    double currentStepSize = integrationSettingsPointer_->getInitialStepSize();

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

    Eigen::MatrixXd dummy_variable;

    // Declarations for variables used in time loop
    Eigen::MatrixXd currentThrustAccelerationVector( 1, 3 );
    double previousTime;
    Eigen::MatrixXd  nextStateVector(  numberOfStateVariables, 1 );
    Eigen::MatrixXd  coefficientMatrix( order + 1, numberOfStateVariables );
    int numberOfRowsInStateMatrix;
    Eigen::MatrixXd interTime( order + 1, 1 );
    Eigen::MatrixXd coefficientMatrixInLoop( order + 1, 1 );


    // Advance the state
    while( currentTime < finalTime ) // while the current time is smaller than the final time
    {


        /// Use discrete force model to update the current force,
        /// force derivatives, acceleration and acceleration derivatives in the spacecraft body
        // Construct discreteForceModel
        DiscreteForceModelPointer discreteForceModelPointer =
                boost::make_shared< DiscreteForceModel >(
                    spacecraftPointer_,
                    integrationSettingsPointer_,
                    constantsPointer_);
        // Update variables in spacecraft body
        discreteForceModelPointer->updateCurrentForcesAndAccelerationsForTSI( currentTime ); // Writes output in spacecraftPointer_
        // is not done for the finalTime. So at finalTime the current forces and accelerations are not calculated


        // Get current thrust accelerationvector from spacecraft pointer. Must be computed before next step size is calculated
        currentThrustAccelerationVector = spacecraftPointer_->getCurrentAcceleration();


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
        problemRecurrenceRelationsPointer_->computeNextState( currentStepSize);

        //! Get output of problemRecurrenceRelations
        // Get current state vector from spacecraft pointer
        nextStateVector = spacecraftPointer_->getCurrentState();
        // Get current coefficients of derivative matrix from spacecraft pointer.
        coefficientMatrix = spacecraftPointer_->getCurrentCoefficientMatrix();

        // Store the next state in new row of the state matrix
        numberOfRowsInStateMatrix = stateMatrix.rows();
        stateMatrix.conservativeResize( numberOfRowsInStateMatrix + 1, Eigen::NoChange );
        stateMatrix.row( numberOfRowsInStateMatrix ) = nextStateVector.adjoint();

        // Store the next time in new row of the time vector
        Time.conservativeResize( numberOfRowsInStateMatrix + 1, Eigen::NoChange );
        Time( numberOfRowsInStateMatrix, 0 ) = currentTime;

        // Store the next computed thrust acceleration vector in a new row of the computed thrust acceleration matrix
        computedThrustAccelerationMatrix.conservativeResize( numberOfRowsInStateMatrix, Eigen::NoChange );
        computedThrustAccelerationMatrix.row( numberOfRowsInStateMatrix - 1 ) = currentThrustAccelerationVector;

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

            // Add refined states to state matrix
            Eigen::MatrixXd temp_stateMatrix( stateMatrix.rows() + refineFactor, stateMatrix.cols() );
            temp_stateMatrix << stateMatrix.topRows( stateMatrix.rows() - 1 ),
                                newInterStateMatrix,
                                stateMatrix.bottomRows( 1 );

            stateMatrix = temp_stateMatrix;

            // Add refined times to state matrix
            Eigen::MatrixXd temp_Time( Time.rows() + refineFactor, Time.cols() );
            temp_Time << Time.topRows( Time.rows() - 1 ),
                         referenceTime,
                         Time.bottomRows( 1 );

            Time = temp_Time;
        }

    }

    // Store output
    setStateMatrix( stateMatrix );
    setTimeVector( Time );
    setComputedThrustAccelerationMatrix( computedThrustAccelerationMatrix );



}


//void TaylorSeriesIntegrator::integrate()
//{
//    // Set loop counter
//    int i = 0; // Starting from the first element in a vector

//    Vector8d errorRelative = settings_.getErrorRelative();
//    Vector8d errorAbsolute = settings_.getErrorAbsolute();

//    Eigen::VectorXd Time;
//    Time( 0 ) = settings_.getInitialTime();

//    // Get thrust force matrix
//    // Put thrustForceMatrix_ in an std vector of std vector
//    std::vector< std::vector > thrustForceMatrix_;


//    std::vector thrustForceMatrix_; // vector of vectors, so it forms a matrix
//    thrustForceMatrix_= spacecraft_.getThrustForce();

//    // Extract thrust force vectors and time vector out of
//    // the thrust force matrix
//    Eigen::VectorXd thrustForceValue_time = thrustForceMatrix_.col(0);
//    Eigen::VectorXd thrustForceVector_e1 = thrustForceMatrix_.col(1);
//    Eigen::VectorXd thrustForceVector_e2 = thrustForceMatrix_.col(2);
//    Eigen::VectorXd thrustForceVector_e3 = thrustForceMatrix_.col(3);

//    while ( Time( i ) < ( settings_.getFinalTime() - settings_.getStepSize() ) )
//    {

//        // Calculation of limitingError
//        Vector8d reducedCurrentState = spacecraft_.getState().block(0,0,8,1); // Extract from element (1,1) to the eighth row and first column
//        Vector8d reducedCurrentStateAbs = reducedCurrentState.cwiseAbs();
//        Vector8d product = errorRelative.cwiseProduct( reducedCurrentStateAbs );
//        Vector8d limitingError = product.cwiseMax( errorAbsolute );


//        //! Interpolator

//        // Extract thrust force vectors and time vector out of
//        // the thrust force matrix has been done outside while loop already

//        const double thrustForceValue_e1 = tudat::interpolators::computeLinearInterpolation(
//                    thrustForceValue_time, thrustForceVector_e1,
//                    Time( i ) );

//        const double thrustForceValue_e2 = tudat::interpolators::computeLinearInterpolation(
//                    thrustForceValue_time, thrustForceVector_e2,
//                    Time( i ) );

//        const double thrustForceValue_e3 = tudat::interpolators::computeLinearInterpolation(
//                    thrustForceValue_time, thrustForceVector_e3,
//                    Time( i ) );

//        // Store the three thrust force values in a 1x3 vector
//        Eigen::Vector3d thrustForceVector_interpolated;
//        thrustForceVector_interpolated( thrustForceValue_e1, thrustForceValue_e2, thrustForceValue_e3 );

//// Shouldn't there also be a second value of the thrust force vector be provided to the integratorStep function??
//        integratorStep( settings_.getStepSize(), thrustForceVector_interpolated, limitingError );

//        ///! Store next time and increase loop counter

//        Time( i + 1 ) = Time( i ) + settings_.getStepSize();
//        i++;
//    }
//}

void TaylorSeriesIntegrator::integratorStep(double& stepSize,const Eigen::Vector3d& thrustForceVector,const double limitingError)
{

}

