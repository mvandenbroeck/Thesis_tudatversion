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

#include "USMPropagation/recurrenceRelations.h"


namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_recurrenceRelations )

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testMultiplyRecursive )
{


    // Define order
    int order = 4;

    // Set order
    RecurrenceRelations testRecurrenceRelation( order ); //

    // Define input
    Eigen::VectorXd X_left( 5 ), U_left( 5 ), X_right( 5 ), U_right( 5 );
    X_left << 0.5, // Must containt order+1 elements
              0.4,
              0.2,
              0.6,
              1.2;
    U_left << 0.2, // Must containt order+1 elements
              1.4,
              6.2,
              -0.6,
              1.3;
    X_right << 0.8, // Must containt order+1 elements
               0.6,
               0.1,
               0.9,
               1.6;
    U_right << 1.5, // Must containt order+1 elements
               2.4,
               1.2,
               -0.6,
               -1.2;

    // Expected output
    double expectedOutput;
    expectedOutput = 3.8250000000; // Calculated with Matlab code

    // Compute output
    double computedOutput;
    computedOutput = testRecurrenceRelation.multiplyRecursive( X_left, U_left, X_right, U_right );

    // Check if computed output matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedOutput, expectedOutput,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testV1overV2 )
{


    // Define order
    int order = 4;

    // Set order
    RecurrenceRelations testRecurrenceRelation( order ); //

    // Define input
    Eigen::VectorXd W( 5 ), V1( 5 ), V2( 5 );
    W << 0.5, // Must containt order+1 elements
              0.4,
              0.2,
              0.6,
              1.2;
    V1 << 0.2, // Must containt order+1 elements
              1.4,
              6.2,
              -0.6,
              1.3;
    V2 << 0.8, // Must containt order+1 elements
               0.6,
               0.1,
               0.9,
               1.6;

    // Expected output
    double expectedOutput;
    expectedOutput = -0.3000000000; // Calculated with Matlab code

    // Compute output
    double computedOutput;
    computedOutput = testRecurrenceRelation.V1overV2( W, V1, V2 );

    // Check if computed output matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedOutput, expectedOutput,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testDivide_Recursive )
{

    // Define order
    int order = 10;

    // Set order
    RecurrenceRelations testRecurrenceRelation( order ); //

    // Define input
    Eigen::VectorXd W( order + 1 ), U_num( order + 1 ),
            U_denom( order + 1 ), X_denom( order + 1 );
    W << 0.5, // Must containt order+1 elements
              0.4,
              0.2,
              0.6,
              1.2,
              4.0,
              6.2,
              -3.82,
              0.01,
              2.3,
              0.88;
    U_num << 0.2, // Must containt order+1 elements
              1.4,
              6.2,
              -0.6,
              1.3,
              0.81,
              -1.2,
              4.25,
              3.01,
              9.0,
              2.5;
    U_denom << 0.8, // Must containt order+1 elements
               0.6,
               0.1,
               0.9,
               1.6,
               0.11,
                -0.2,
                -4.25,
                -3.01,
                -9.0,
                4.33;
    X_denom << 1.8, // Must containt order+1 elements
               -0.6,
               1.1,
               -0.9,
               2.6,
               -0.11,
                1.2,
                4.25,
                2.01,
                9.0,
                7.5;

    // Expected output
    double expectedOutput;
    expectedOutput = -1.558608906525573; // Calculated with Matlab code

    // Compute output
    double computedOutput;
    computedOutput = testRecurrenceRelation.divideRecursive( W, U_num, U_denom, X_denom );

    // Check if computed output matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedOutput, expectedOutput,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testUm_over_Xn )
{

    // Define order
    int order = 10;

    // Set order
    RecurrenceRelations testRecurrenceRelation( order ); //

    // Define input
    Eigen::VectorXd W( order + 1 ), U_num( order + 1 ),
            U_denom( order + 1 ), X_denom( order + 1 );
    W << 0.5, // Must containt order+1 elements
              0.4,
              0.2,
              0.6,
              1.2,
              4.0,
              6.2,
              -3.82,
              0.01,
              2.3,
              0.01;
    U_num << 0.2, // Must containt order+1 elements
              1.4,
              6.2,
              -0.6,
              1.3,
              0.81,
              -1.2,
              4.25,
              3.01,
              9.0,
              5.2;
    U_denom << 0.8, // Must containt order+1 elements
               0.6,
               0.1,
               0.9,
               1.6,
               0.11,
                -0.2,
                -4.25,
                -3.01,
                -9.0,
                0.6;
    X_denom << 1.8, // Must containt order+1 elements
               -0.6,
               1.1,
               -0.9,
               2.6,
               -0.11,
                1.2,
                4.25,
                2.01,
                9.0,
                1.1;

    // Expected output
    double expectedOutput;
    expectedOutput = 0.83027998236331584; // Calculated with Matlab code

    // Compute output
    double computedOutput;
    computedOutput = testRecurrenceRelation.Um_over_Xn( W, U_num, U_denom, X_denom );

    // Check if computed output matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedOutput, expectedOutput,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testX1xV2_Recursive )
{

    // Define order
    int order = 10;

    // Set order
    RecurrenceRelations testRecurrenceRelation( order ); //

    // Define input
    Eigen::VectorXd X_L( order + 1 ), U_L( order + 1 ),
            V_R( order + 1 );
    X_L << 0.5, // Must containt order+1 elements
              0.4,
              0.2,
              0.6,
              1.2,
              4.0,
              6.2,
              -3.82,
              0.01,
              2.3,
              0.01;
    U_L << 0.2, // Must containt order+1 elements
              1.4,
              6.2,
              -0.6,
              1.3,
              0.81,
              -1.2,
              4.25,
              3.01,
              9.0,
              5.2;
    V_R << 0.8, // Must containt order+1 elements
               0.6,
               0.1,
               0.9,
               1.6,
               0.11,
                -0.2,
                -4.25,
                -3.01,
                -9.0,
                0.6;

    // Expected output
    double expectedOutput;
    expectedOutput =  -11.296227380952383; // Calculated with Matlab code

    // Compute output
    double computedOutput;
    computedOutput = testRecurrenceRelation.X1xV2_Recursive( X_L, U_L, V_R );

    // Check if computed output matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedOutput, expectedOutput,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testV1xV2 )
{

    // Define order
    int order = 10;

    // Set order
    RecurrenceRelations testRecurrenceRelation( order ); //

    // Define input
    Eigen::VectorXd V1( order + 1 ), V2( order + 1 );
    V1 << 0.5, // Must containt order+1 elements
              0.4,
              0.2,
              0.6,
              1.2,
              4.0,
              6.2,
              -3.82,
              0.01,
              2.3,
              0.01;
    V2 << 0.2, // Must containt order+1 elements
              1.4,
              6.2,
              -0.6,
              1.3,
              0.81,
              -1.2,
              4.25,
              3.01,
              9.0,
              5.2;

    // Expected output
    double expectedOutput;
    expectedOutput =  24.788000000000000; // Calculated with Matlab code

    // Compute output
    double computedOutput;
    computedOutput = testRecurrenceRelation.V1xV2( V1, V2 );

    // Check if computed output matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedOutput, expectedOutput,
                                std::numeric_limits< double >::epsilon( ) );
}


//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testX1xU2overX2_Recursive )
{

    // Define order
    int order = 10;

    // Set order
    RecurrenceRelations testRecurrenceRelation( order ); //

    // Define input
    Eigen::VectorXd W( order + 1 ), X_L( order + 1 ),
            U_L( order + 1 ), X_R( order + 1 ), U_R( order + 1 );
    W << 0.5, // Must containt order+1 elements
              0.4,
              0.2,
              0.6,
              1.2,
              4.0,
              6.2,
              -3.82,
              0.01,
              2.3,
              0.01;
    X_L << 0.2, // Must containt order+1 elements
              1.4,
              6.2,
              -0.6,
              1.3,
              0.81,
              -1.2,
              4.25,
              3.01,
              9.0,
              5.2;
    U_L << 0.8, // Must containt order+1 elements
               0.6,
               0.1,
               0.9,
               1.6,
               0.11,
                -0.2,
                -4.25,
                -3.01,
                -9.0,
                0.6;
    X_R << 1.8, // Must containt order+1 elements
               -0.6,
               1.1,
               -0.9,
               2.6,
               -0.11,
                1.2,
                4.25,
                2.01,
                9.0,
                1.1;
    U_R << -3.8, // Must containt order+1 elements
               2.6,
               -2.1,
               2.9,
               -4.6,
               2.11,
                -3.2,
                -6.25,
                -4.01,
                -11.0,
                -3.1;


    // Expected output
    double expectedOutput;
    expectedOutput = -0.7791646825396844633; // Calculated with Matlab code

    // Compute output
    double computedOutput;
    computedOutput = testRecurrenceRelation.X1xU2overX2_Recursive( W, X_L, U_L, X_R, U_R );

    // Check if computed output matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedOutput, expectedOutput,
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
