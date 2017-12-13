#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Core>
#include "spacecraft.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "constants.h"

using namespace std;

Spacecraft::Spacecraft(
        double initialMass,
        Vector7d initialReducedState,
        double specificImpulse,
        ConstantsPointer constantsPointer ) :
    mass_( initialMass ),
    specificImpulse_( specificImpulse ),
    constantsPointer_( constantsPointer ),
    initialReducedState_( initialReducedState )
{
    for ( int i = 0; i < initialReducedState_.size(); i++ ) // dim( initialReducedState ) = 7x1
    {
        initialState_( i ) = initialReducedState_( i );
    }
    double cosineLambda = ( pow(initialState_( 6 ),2) - pow(initialState_( 5 ),2) )
            / ( pow(initialState_( 5 ),2) + pow(initialState_( 6 ),2) );
    double sineLambda = ( 2 * initialState_( 5 ) * initialState_( 6 ) )
            / ( pow(initialState_( 5 ),2) + pow(initialState_( 6 ),2) );
    // mass is already included in initialReducedState vector
    initialState_( 7 ) = mass_;
    initialState_( 8 ) = initialState_( 1 ) * cosineLambda + initialState_( 2 ) * sineLambda;
    initialState_( 9 ) = initialState_( 0 ) - initialState_( 1 ) * sineLambda
            + initialState_( 2 ) * cosineLambda;
    initialState_( 10 ) = initialState_( 0 ) * pow(initialState_( 9 ),2)
            / constantsPointer_->centralBodyGravitationalParameter_;
    currentState_ = initialState_;
}


//Spacecraft::Spacecraft(array<double,8>& reducedState, double specificImpulse ) : specificImpulse_(specificImpulse)
//{

//    for (size_t i = 0; i < reducedState.size(); i++)
//    {
//        initialState_[i] = reducedState[i];
//    }
//    double cosineLambda = ( pow(initialState_[6],2) - pow(initialState_[5],2) ) / ( pow(initialState_[5],2) + pow(initialState_[6],2) );
//    double sineLambda = ( 2 * initialState_[5] * initialState_[6] ) / ( pow(initialState_[5],2) + pow(initialState_[6],2) );
//    initialState_[8] = initialState_[1] * cosineLambda + initialState_[2] * sineLambda;
//    initialState_[9] = initialState_[0] - initialState_[1] * sineLambda + initialState_[2] * cosineLambda;
//    initialState_[10] = initialState_[0] * pow(initialState_[9],2) / constants_.getCentralBodyGravitationalParameter();
//    state_ = initialState_;
//}
//void Spacecraft::loadThrustForce(string fileName)
//{
//    fstream thrustForceFile (fileName,ios::in);
//    if (!thrustForceFile)
//    {
//        cerr << "File could not be opened" << endl;
//        exit(EXIT_FAILURE);
//    }
//    int row = 0;
//    array<double,4> dummy;
//    thrustForce_.push_back(dummy);
//    while (thrustForceFile >> thrustForce_[row][0] >> thrustForce_[row][1]
//           >> thrustForce_[row][2] >> thrustForce_[row][3])
//    {
//        thrustForce_.push_back(dummy);
//        row++;
//    }
//    thrustForce_.pop_back();
//}

void Spacecraft::loadThrustForce( string fileName)
{
    // Read input file and store data in matrix.
    thrustForceMatrix_ = tudat::input_output::readMatrixFromFile( fileName , " \t", "#" );
}
