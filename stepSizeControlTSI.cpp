#include <cmath>
#include <iostream>
#include "stepSizeControlTSI.h"
#include "integrationSettings.h"

StepSizeControlTSI::StepSizeControlTSI( IntegrationSettingsPointer integrationSettingsPointer ):
    integrationSettingsPointer_( integrationSettingsPointer ) {}

double StepSizeControlTSI::determineNextStepSize(
        Eigen::MatrixXd derivatives, double previousStepSize )
{
    // Get the order index and the amount of state variables
    int K, n_State;
    K = derivatives.rows() - 1;
    n_State = derivatives.cols();

    Eigen::MatrixXd X( K + 1, 1);

    Eigen::MatrixXd h_next_vect( n_State, 1 );

    double nextStepSize = previousStepSize;

    do
    {
        previousStepSize = nextStepSize;

        for( int n = 0; n < n_State; n = n + 1 )
        {
            X = derivatives.col(n);
            h_next_vect( n ) =  std::exp( ( 1.0 / ( K - 1 ) ) *
                                          std::log(  integrationSettingsPointer_->getErrorTolerance( n )
                                                     / ( std::abs( X( K - 1 ) )
                                                         + previousStepSize * std::abs( X( K ) ) ) ) );
        }
        // Next step-size is smallest demanded stepsize
        nextStepSize = h_next_vect.minCoeff();


        //??
//        std::cout << "previousStepSize = " << previousStepSize << std::endl;
//        std::cout << "nextStepSize = " << nextStepSize << std::endl;
//        std::cout << "Loop condition = " << ( std::abs( nextStepSize - previousStepSize ) / previousStepSize ) << std::endl;

    } while ( ( std::abs( nextStepSize - previousStepSize ) / previousStepSize ) > 1.0e-6 );

    nextStepSize = integrationSettingsPointer_->getSafetyFactor() * nextStepSize;

    //??
//    std::cout << "Step with size " << nextStepSize << std::endl;

    return nextStepSize;


}

StepSizeControlTSI::~StepSizeControlTSI()
{

}

