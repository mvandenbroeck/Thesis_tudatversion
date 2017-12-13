#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Core>
#include "spacecraftGen.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "constants.h"

using namespace std;

SpacecraftGen::SpacecraftGen(
        double initialMass,
        Eigen::VectorXd initialReducedState,
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
    // mass is already included in initialReducedState vector
    initialState_( 7 ) = mass_;
    currentState_ = initialState_;
}

void SpacecraftGen::loadThrustForce( string fileName)
{
    // Read input file and store data in matrix.
    thrustForceMatrix_ = tudat::input_output::readMatrixFromFile( fileName , " \t", "#" );
}
