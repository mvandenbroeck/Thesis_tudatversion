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
 *      251017    M. Van den Broeck Creation.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <Eigen/Core>

//#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Thesis/lagrangeInterpolator.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_lagrangeInterpolator )

//! Test if next step size is computed correctly.
BOOST_AUTO_TEST_CASE( test1 )
{
    // Define points to be interpolated
    Eigen::Matrix< double, 3, 2 > pointsToBeInterpolated;
    pointsToBeInterpolated << -1.0, -2.0,
                               0.0, -1.0,
                               1.0,  0.0;


    // Define expected output
    Eigen::Matrix< double, 3, 1> expectedPolynomialCoefficients;
    expectedPolynomialCoefficients << 0.0,
                                      1.0,
                                     -1.0;

    // Construct lagrange interpolator using pointer
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer =
            boost::make_shared< LagrangeInterpolator >( pointsToBeInterpolated );

    // Determine next step size
    Eigen::Matrix< double, 3, 1> computedPolynomialCoefficients;
    computedPolynomialCoefficients = lagrangeInterpolatorPointer->computeCoefficientsOfInterpolatingPolynomial();

    // Check if computed next step size matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 0, 0 ), expectedPolynomialCoefficients( 0, 0),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 1, 0 ), expectedPolynomialCoefficients( 1, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 2, 0 ), expectedPolynomialCoefficients( 2, 0 ),
                                1.0e-15 );
    //std::numeric_limits< double >::epsilon( )

}

//! Test if the coefficients of the 2nd order Lagrange Polynomial are computed correctly.
BOOST_AUTO_TEST_CASE( test2 )
{
    // Define points to be interpolated
    Eigen::Matrix< double, 3, 2 > pointsToBeInterpolated;
    pointsToBeInterpolated << -10.0, -2.0,
                               0.0, -1.0,
                               1.0,  10.0;


    // Define expected output
    Eigen::Matrix< double, 3, 1> expectedPolynomialCoefficients;
    expectedPolynomialCoefficients << 0.9909090909090909,
                                      10.009090909090909,
                                     -1.0;

    // Construct lagrange interpolator using pointer
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer =
            boost::make_shared< LagrangeInterpolator >( pointsToBeInterpolated );

    // Determine next step size
    Eigen::Matrix< double, 3, 1> computedPolynomialCoefficients;
    computedPolynomialCoefficients = lagrangeInterpolatorPointer->computeCoefficientsOfInterpolatingPolynomial();

    // Check if computed next step size matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 0, 0 ), expectedPolynomialCoefficients( 0, 0),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 1, 0 ), expectedPolynomialCoefficients( 1, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 2, 0 ), expectedPolynomialCoefficients( 2, 0 ),
                                1.0e-15 );
    //std::numeric_limits< double >::epsilon( )

}

//! Test if the coefficients of the 2nd order Lagrange Polynomial are computed correctly.
BOOST_AUTO_TEST_CASE( test3 )
{
    // Define points to be interpolated
    Eigen::Matrix< double, 3, 2 > pointsToBeInterpolated;
    pointsToBeInterpolated << -10.0, -2.0,
                               0.0, 1.0,
                               5.0,  10.0;


    // Define expected output
    Eigen::Matrix< double, 3, 1> expectedPolynomialCoefficients;
    expectedPolynomialCoefficients << 0.1,
                                      1.3,
                                      1.0;

    // Construct lagrange interpolator using pointer
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer =
            boost::make_shared< LagrangeInterpolator >( pointsToBeInterpolated );

    // Determine next step size
    Eigen::Matrix< double, 3, 1> computedPolynomialCoefficients;
    computedPolynomialCoefficients = lagrangeInterpolatorPointer->computeCoefficientsOfInterpolatingPolynomial();

    // Check if computed next step size matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 0, 0 ), expectedPolynomialCoefficients( 0, 0),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 1, 0 ), expectedPolynomialCoefficients( 1, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 2, 0 ), expectedPolynomialCoefficients( 2, 0 ),
                                1.0e-15 );
    //std::numeric_limits< double >::epsilon( )

}

//! Test if the coefficients of the 3rd order Lagrange Polynomial are computed correctly.
BOOST_AUTO_TEST_CASE( test4 )
{
    // Define points to be interpolated
    Eigen::Matrix< double, 4, 2 > pointsToBeInterpolated;
    pointsToBeInterpolated << -10.0,  -2.0,
                                0.0,   1.0,
                                5.0,  10.0,
                                6.0,  12.0;


    // Define expected output
    Eigen::Matrix< double, 4, 1> expectedPolynomialCoefficients;
    expectedPolynomialCoefficients << -0.0041666666666666666,
                                      0.079166666666666666,
                                     1.50833333333333333,
                                     1.0;

    // Construct lagrange interpolator using pointer
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer =
            boost::make_shared< LagrangeInterpolator >( pointsToBeInterpolated );

    // Determine next step size
    Eigen::Matrix< double, 4, 1> computedPolynomialCoefficients;
    computedPolynomialCoefficients = lagrangeInterpolatorPointer->computeCoefficientsOfInterpolatingPolynomial();

    // Check if computed next step size matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 0, 0 ), expectedPolynomialCoefficients( 0, 0),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 1, 0 ), expectedPolynomialCoefficients( 1, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 2, 0 ), expectedPolynomialCoefficients( 2, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 3, 0 ), expectedPolynomialCoefficients( 3, 0 ),
                                1.0e-15 );
    //std::numeric_limits< double >::epsilon( )

}

//! Test if the coefficients of the 3rd order Lagrange Polynomial are computed correctly.
BOOST_AUTO_TEST_CASE( test5 )
{
    // Define points to be interpolated
    Eigen::Matrix< double, 4, 2 > pointsToBeInterpolated;
    pointsToBeInterpolated << 0.0,  5.0,
                              86400.0,   6.0,
                              2.0 * 86400.0,  2.0,
                              3.0 * 86400.0,  0.0;


    // Define expected output
    Eigen::Matrix< double, 4, 1> expectedPolynomialCoefficients;
    expectedPolynomialCoefficients << 1.80886252836627e-15,
                                      -8.03755144032922e-10,
                                     6.75154320987654e-05,
                                     5;

    // Construct lagrange interpolator using pointer
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer =
            boost::make_shared< LagrangeInterpolator >( pointsToBeInterpolated );

    // Determine next step size
    Eigen::Matrix< double, 4, 1> computedPolynomialCoefficients;
    computedPolynomialCoefficients = lagrangeInterpolatorPointer->computeCoefficientsOfInterpolatingPolynomial();

    // Check if computed next step size matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 0, 0 ), expectedPolynomialCoefficients( 0, 0),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 1, 0 ), expectedPolynomialCoefficients( 1, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 2, 0 ), expectedPolynomialCoefficients( 2, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedPolynomialCoefficients( 3, 0 ), expectedPolynomialCoefficients( 3, 0 ),
                                1.0e-15 );
    //std::numeric_limits< double >::epsilon( )

}

//! Test if the coefficients of the 3rd order Lagrange Polynomial are computed correctly.
BOOST_AUTO_TEST_CASE( test6 )
{
    // Define points to be interpolated
    Eigen::Matrix< double, 4, 2 > pointsToBeInterpolated;
    pointsToBeInterpolated << -10.0,  -2.0,
                                0.0,   1.0,
                                5.0,  10.0,
                                6.0,  12.0;


    // Define expected output
    Eigen::Matrix< double, 4, 1> expectedValue;
    expectedValue = pointsToBeInterpolated.col( 1 );

    // Construct lagrange interpolator using pointer
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer =
            boost::make_shared< LagrangeInterpolator >( pointsToBeInterpolated );

    // Determine next step size
    Eigen::Matrix< double, 4, 1> computedPolynomialCoefficients;
    computedPolynomialCoefficients = lagrangeInterpolatorPointer->computeCoefficientsOfInterpolatingPolynomial();

    double currentTime0 = -10.0;

    Eigen::Matrix< double, 4, 1> t_powers0;
    t_powers0 <<   currentTime0 * currentTime0 * currentTime0,
                   currentTime0 * currentTime0,
                   currentTime0,
                   1.0;

    double currentTime1 = 0.0;

    Eigen::Matrix< double, 4, 1> t_powers1;
    t_powers1 <<   currentTime1 * currentTime1 * currentTime1,
                   currentTime1 * currentTime1,
                   currentTime1,
                   1.0;

    double currentTime2 = 5.0;

    Eigen::Matrix< double, 4, 1> t_powers2;
    t_powers2 <<   currentTime2 * currentTime2 * currentTime2,
                   currentTime2 * currentTime2,
                   currentTime2,
                   1.0;

    double currentTime3 = 6.0;

    Eigen::Matrix< double, 4, 1> t_powers3;
    t_powers3 <<   currentTime3 * currentTime3 * currentTime3,
                   currentTime3 * currentTime3,
                   currentTime3,
                   1.0;

    // Computed interpolated values at nodes
    double computedValue0, computedValue1, computedValue2, computedValue3;
    computedValue0 = computedPolynomialCoefficients.adjoint() * t_powers0;
    computedValue1 = computedPolynomialCoefficients.adjoint() * t_powers1;
    computedValue2 = computedPolynomialCoefficients.adjoint() * t_powers2;
    computedValue3 = computedPolynomialCoefficients.adjoint() * t_powers3;

    // Check if computed interpolated values match the expected values.
    BOOST_CHECK_CLOSE_FRACTION( computedValue0, expectedValue( 0, 0),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedValue1, expectedValue( 1, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedValue2, expectedValue( 2, 0 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedValue3, expectedValue( 3, 0 ),
                                1.0e-15 );
    //std::numeric_limits< double >::epsilon( )

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
