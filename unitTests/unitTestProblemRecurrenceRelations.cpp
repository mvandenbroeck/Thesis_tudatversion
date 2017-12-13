/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      161017    M. Van den Broeck Creation.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <Eigen/Core>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Thesis/constants.h"
#include "Thesis/spacecraft.h"
#include "Thesis/recurrenceRelations.h"
#include "Thesis/problemRecurrenceRelations.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_problemRecurrenceRelations )

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( test1 )
{

    // Define order
    int order = 4;

    // Define constants
    double standardGravity = 9.81;
    double centralBodyGravitationalParameter = 3.986004418e14; // Test case = orbit around Earth

    // // Create constants object
    ConstantsPointer constantsPointer =
            boost::make_shared< Constants >( centralBodyGravitationalParameter, standardGravity  );

    // // Create spacecraft object
    double specificImpulse = 3000.0;
    // USM state
    Vector7d initialReducedState;
    initialReducedState << 7433.80084166613,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0;
    double initialMassOfSpacecraft = 2000.0;

    SpacecraftPointer spacecraftPointer =
            boost::make_shared< Spacecraft >( initialMassOfSpacecraft,
                                              initialReducedState,
                                              specificImpulse,
                                              constantsPointer );

    // Set current acceleration and its derivatives
    Eigen::Matrix< double, 1, 3 > acc_curr;
    acc_curr << 0.0, 0.005, 0.0;
    spacecraftPointer->setCurrentAcceleration( acc_curr );

    Eigen::MatrixXd dacc_curr( order, 3 );
    dacc_curr.fill( 0.0 );
    dacc_curr( 0, 1 ) = 4.12422538201892e-07; // From Matlab code
    dacc_curr( 1, 1 ) = -1.3565332505994e-10; // From Matlab code
    spacecraftPointer->setCurrentAccelerationDerivative( dacc_curr );

    // Set derivatives of current resultant thrust force.
    Eigen::MatrixXd dT_res( order, 1 );
    dT_res << 0.000823994726433751,
              -2.71586923676253e-07,
              0.0,
              0.0;
    spacecraftPointer->setCurrentResultantThrustForceDerivatives( dT_res );


    // // Create object using constructor

    ProblemRecurrenceRelationsPointer problemRecurrenceRelationPointer =
            boost::make_shared< ProblemRecurrenceRelations>( order, spacecraftPointer, constantsPointer );


    // // Define input for computeState function

    // Step size
    double currentStepSize = 1.0;

    // Next state vector
    Eigen::MatrixXd nextStateVector( 11, 1); // Will be overwritten inside function

    // Next coefficient matrix
    Eigen::MatrixXd coefficientMatrix( 5, 11 ); // Will be overwritten inside function


    // Set order
    problemRecurrenceRelationPointer->computeNextState( currentStepSize ); //

    // Expected output
    Eigen::Matrix< double, 11, 1 > expectedNextStateVector;
    expectedNextStateVector << 7433.79584146333, // Calculated with Matlab code
                               -5.15333924737626e-06,
                               0.0100004071984344,
                               0.0,
                               0.0,
                               0.000515305906886064,
                               0.999999867229902,
                               1999.99966019667,
                               5.15319992233846e-06,
                               7433.80584187053,
                               0.00103061220599872;

    Eigen::Matrix< double, 5, 11 > expectedCoefficientMatrix; // Calculated with Matlab code                        //                                                                          //                                                         //
    expectedCoefficientMatrix << 7433.80084166613,                  0.0,                    0.0,                    0.0,                        0.0,                    0.0,                    1.0,                    2000.0,                 0.0,       7433.80084166613,       0.00103061151277778,
                                 -0.005,                            0.0,                   0.01,                    0.0,                        0.0,    0.00051530575638889,                    0.0,     -0.000339789330615019,                 0.0,                  0.005,      6.93192846249827e-10,
                                 -2.02848251913308e-07, -5.1530575638889e-06, 4.09059521014254e-07,                 0.0,                        0.0,   1.73298211562457e-10,  -1.32770011283763e-07,     -1.39992308262615e-08, 5.1530575638889e-06,   2.06211269100946e-07,      2.85888353113748e-14,
                                 4.54929098599676e-11, -2.8220962258952e-10, -1.86097750199679e-09,                 0.0,                        0.0,  -2.28009522242284e-11,  -8.93015659900336e-14,       1.5380389833291e-12,1.42837793410384e-10,  -1.81548459213682e-09,     -4.97123445064975e-16,
                                 -3.53057259164559e-16, 5.26135234416357e-13, -1.09161632842322e-13,                0.0,                        0.0,  -8.51492361377252e-17,   2.93550896801263e-15,                       0.0,-4.79343848272888e-13, -3.76955842470294e-14,     -1.04933925385682e-20;



    // Computed output
    Eigen::Matrix< double, 11, 1 > computedNextStateVector;
    computedNextStateVector = spacecraftPointer->getCurrentState();

    Eigen::Matrix< double, 5, 11 > computedCoefficientMatrix;
    computedCoefficientMatrix = spacecraftPointer->getCurrentCoefficientMatrix();

    // Check if computed sample mean matches expected value.
    for ( int i = 0; i < computedNextStateVector.rows(); i++ )
    {
        //std::cout << "i =  " << i << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedNextStateVector( i ), expectedNextStateVector( i ),
                                4.9e-15 );
    }

    for ( int i = 0; i < ( computedCoefficientMatrix.rows() * computedCoefficientMatrix.cols() ); i++ )
    {
        //std::cout << "i = " << i <<  std::endl;
        BOOST_CHECK_CLOSE_FRACTION( computedCoefficientMatrix( i ), expectedCoefficientMatrix( i ),
                                4.9e-15 );
    }

}




BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
