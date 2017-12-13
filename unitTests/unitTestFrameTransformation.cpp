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
#include <iostream>

//#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Thesis/frameTransformation.h"

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_frameTransformation )

//! Test if next step size is computed correctly.
BOOST_AUTO_TEST_CASE( test1 )
{
    //??
    std::cout << "Label 0" << std::endl;

    // Define vehicle state in Cartesian elements
    Eigen::Matrix< double, 6, 1 > vehicleStateCartesian;
    //??
    std::cout << "Label 1" << std::endl;
    vehicleStateCartesian( 0, 0 ) = 0.0;
    vehicleStateCartesian( 1, 0 ) = 4.0;
    vehicleStateCartesian( 2, 0 ) = 4.0;
    vehicleStateCartesian( 3, 0 ) = 15.0;
    vehicleStateCartesian( 4, 0 ) = 0.0;
    vehicleStateCartesian( 5, 0 ) = 0.0;

    //??
    std::cout << "Label 2" << std::endl;

    // Define thrust in velocity frame
    Eigen::Matrix< double, 1, 3 > thrustInVelocityFrame;
    thrustInVelocityFrame <<
            0.0,
            2.0,
            0.0;

    // Define expected output
    Eigen::Matrix< double, 1, 3 > expectedThrustInCartFrame;
    expectedThrustInCartFrame <<
            0.0,
            -std::sqrt(2.0),
            -std::sqrt(2.0);

    // Construct frameTransformation object
    FrameTransformation frameTransformation;

    // Compute thrust in Cartesian frame
    Eigen::Matrix< double, 1, 3 > computedThrustInCartFrame;
    computedThrustInCartFrame = frameTransformation.velocityFrameToCartFrame(
                thrustInVelocityFrame,
                vehicleStateCartesian );

    std::cout << "computedThrustInCartFrame =\n" << computedThrustInCartFrame << std::endl;


    // Check if computed thrust in Cartesian frame matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedThrustInCartFrame( 0 ), expectedThrustInCartFrame( 0),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedThrustInCartFrame( 1 ), expectedThrustInCartFrame( 1 ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( computedThrustInCartFrame( 2 ), expectedThrustInCartFrame( 2 ),
                                1.0e-15 );
    //std::numeric_limits< double >::epsilon( )

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
