#include "discreteForceModelGen.h"

#include <Eigen/Core>

#include <boost/make_shared.hpp>

#include "spacecraft.h"
#include "integrationSettings.h"
#include "constants.h"
#include "lagrangeInterpolator.h"
#include "frameTransformation.h"
#include "tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"
#include "tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

//??
#include <iostream>

DiscreteForceModel::DiscreteForceModel(
        SpacecraftPointer spacecraftPointer,
        IntegrationSettingsPointer integrationSettingsPointer,
        ConstantsPointer constantsPointer ):
        spacecraftPointer_( spacecraftPointer ),
        integrationSettingsPointer_( integrationSettingsPointer ),
        constantsPointer_( constantsPointer ) {}

DiscreteForceModel::~DiscreteForceModel()
{

}

void DiscreteForceModelL2::updateCurrentForcesAndAccelerationsForTSI(
        double currentTime)
{
    /// Extract variables from classes
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix;
    thrustForceMatrix = spacecraftPointer_->getThrustForceMatrix();

    /// Specific relations for the specific low-thrust representation method
    // Find index in thrust profile that corresponds with the current time
    int index = tudat::basic_mathematics::computeNearestNeighborUsingBinarySearch(
                thrustForceMatrix.col( 0 ), currentTime);

    // Declarations
    double Time_prev;
    double Time_0;
    double Time_next;

    // First calculate previous, 0 and next times in thrust profile to
    // simplify the next series of calculations
    Time_0 = thrustForceMatrix( index, 0 );
    // When index = 0 -> exception, use first point twice
    if ( index == 0 )
    {
        // The difference between the first and second points is subtracted from first point and the result used as the new previous point
        Time_prev = Time_0 - ( thrustForceMatrix( index + 1, 0 ) - Time_0 );
    }
    else
    {
        Time_prev = thrustForceMatrix( index - 1, 0 ); // scalar of resulting thrust force
    }

    // When index = length( T ) -> exception, use last point twice
    if ( index == ( thrustForceMatrix.rows() - 1 ) )
    {
        // The difference between the last and second-to-last point is added to the last point and the result is used as the new next point
        Time_next = Time_0 + ( Time_0 - thrustForceMatrix( index - 1, 0 ) );
    }
    else
    {
        Time_next = thrustForceMatrix( index + 1, 0 );
    }




    /////// Calculation of derivatives of resulting thrust at current time

    // Calculate derivative of the resulting current thrust force with the use
    // coefficients of the Lagrange polynomial

    // Declarations
    double T_res_prev;
    double T_res_0;
    double T_res_next;

    T_res_0 = ( ( thrustForceMatrix.rightCols( 3 ) ).row( index ) ).norm();
    // When index == 0 -> exception, use first point twice
    if ( index == 0 )
    {
        T_res_prev = T_res_0; // The value of the first point is used twice (assumption)
    }
    else
    {
        T_res_prev = thrustForceMatrix.rightCols( 3 ).row( index - 1 ).norm(); // scalar of resulting thrust force
    }
    // When index = length( T ) -> exception, use last point twice
    if ( index == ( thrustForceMatrix.rows() - 1 )  )
    {
        T_res_next = T_res_0; // The value of the last point is used twice ( assumption )
    }
    else
    {
        T_res_next = thrustForceMatrix.rightCols( 3 ).row( index + 1 ).norm();
    }

    // Calculate the coefficient of the Lagrange polynomial that goes through the
    // points of resulting thrust force
    Eigen::Matrix< double, 3, 2 > PointsToBeInterpolated;
    PointsToBeInterpolated << Time_prev, T_res_prev,
                              Time_0,    T_res_0,
                              Time_next, T_res_next;

    // Construct Lagrange Interpolator
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_T_res =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated );
    Eigen::Matrix< double, 3, 1> lagrangeCoefficients_T_res;
    lagrangeCoefficients_T_res = lagrangeInterpolatorPointer_T_res->computeCoefficientsOfInterpolatingPolynomial();

    // Calculate k-th derivative of the resulting thrust force at the current
    // time. The outcome will be used further on to calculate U8.
    int order = integrationSettingsPointer_->getOrderOfTaylorSeries();
    Eigen::MatrixXd dT_res( std::max( order, 3 ), 1 ); // At least 3 derivatives are calculated
    dT_res.fill( 0.0 );
    // Declarations
    Eigen::Matrix< double, 3, 1> dt_powers_1;
    Eigen::Matrix< double, 3, 1> dt_powers_2;
    // First derivative
    dt_powers_1 << 2.0 * currentTime,
                   1.0,
                   0.0;
    dT_res( 0, 0 ) = lagrangeCoefficients_T_res.adjoint() * dt_powers_1; // First derivative of resulting thrust force
    // Second derivative
    dt_powers_2 << 2.0,
                   0.0,
                   0.0;
    dT_res( 1, 0 ) = lagrangeCoefficients_T_res.adjoint() * dt_powers_2; // Second derivative of resulting thrust force
    // Third to K-th derivatives are zero

    // Store derivatives of resultant low-thrust force in spacecraft class
    spacecraftPointer_->setCurrentResultantThrustForceDerivatives( dT_res );


    /// Calculation of derivatives of acceleration at current time

    // _curr means 'at current step of the integrator'
    // _0 means 'at nearest step in the thrust profile'

    // Extract mass from current state
    Eigen::MatrixXd currentState = spacecraftPointer_->getCurrentState();
    double currentMass = currentState( 7 );

    // Calculate mass flow at nearest step in the thrust profile.
    // Assume that massDerivative_0 = currentMassDerivative (assumption)
    double massDerivative_0;
    massDerivative_0 = - T_res_0 / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );

    // Declarations
    double mass_prev, mass_0, mass_next;
    // Calculate mass at previous and next steps of thrust profile
    mass_prev = currentMass + massDerivative_0 * ( Time_prev - currentTime );
    mass_0 = currentMass + massDerivative_0 * ( Time_0 - currentTime );
    mass_next = currentMass + massDerivative_0 * ( Time_next - currentTime );

    // Transform thrust in velocity frame T_vf to thrust in USM frame
    double velocity_e1 = currentState( 8 );
    double velocity_e2 = currentState( 9 );

    // Calculate flight path angle gamma
    double sin_gamma, cos_gamma, gamma;
    sin_gamma = velocity_e1 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    cos_gamma = velocity_e2 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    gamma = std::atan2( sin_gamma, cos_gamma );

    // Calculate thrust vector at three consecutive points
    Eigen::Matrix< double, 1, 3 > thrust_USM_prev, thrust_USM_0, thrust_USM_next;
    // Create frameTransformation object
    FrameTransformation frameTransformation;
    thrust_USM_0 = frameTransformation.velocityFrameToUSMFrame( thrustForceMatrix.rightCols( 3 ).row( index ), gamma );
    // When index == 0 -> exception, use zero value
    if ( index == 0 )
    {
        thrust_USM_prev = thrust_USM_0; // The value of the first point is used twice (assumption)
    }
    else
    {
        thrust_USM_prev = frameTransformation.velocityFrameToUSMFrame( thrustForceMatrix.rightCols( 3 ).row( index - 1 ), gamma );
    }
    // When index == length( T ) -> exception, use zero value
    if ( index == ( thrustForceMatrix.rows() - 1 ) )
    {
        thrust_USM_next = thrust_USM_0; // The value of the last point is used twice (assumption)
    }
    else
    {
        thrust_USM_next = frameTransformation.velocityFrameToUSMFrame( thrustForceMatrix.rightCols( 3 ).row( index + 1 ), gamma );
    }

    // Calculation of the acceleration vectors at the previous, current and next steps of the thrust profile
    Eigen::Matrix< double, 3, 1 > acc_prev, acc_0, acc_next;
    acc_0 = thrust_USM_0 / mass_0; // 3x1 vector
    acc_prev = thrust_USM_prev / mass_prev; // 3x1 vector
    acc_next = thrust_USM_next / mass_next; // 3x1 vector

    // Calculate the coefficients of the Lagrange polynomial that goes
    // through the three points of the acceleration vector for each
    // direction (e1, e2 and e3).
    Eigen::Matrix< double, 3, 2 > PointsToBeInterpolated_acc_e1,
                                  PointsToBeInterpolated_acc_e2,
                                  PointsToBeInterpolated_acc_e3;
    PointsToBeInterpolated_acc_e1 << Time_prev, acc_prev( 0),
                                     Time_0,    acc_0( 0 ),
                                     Time_next, acc_next( 0 );
    PointsToBeInterpolated_acc_e2 << Time_prev, acc_prev( 1 ),
                                     Time_0,    acc_0( 1 ),
                                     Time_next, acc_next( 1 );
    PointsToBeInterpolated_acc_e3 << Time_prev, acc_prev( 2 ),
                                     Time_0,    acc_0( 2 ),
                                     Time_next, acc_next( 2 );

    // Lagrange interpolation to calculate coefficients
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_acc_e1 =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated_acc_e1 );
    Eigen::Matrix< double, 3, 1> lagrangeCoefficients_acc_e1;
    lagrangeCoefficients_acc_e1 = lagrangeInterpolatorPointer_acc_e1->computeCoefficientsOfInterpolatingPolynomial();

    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_acc_e2 =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated_acc_e2 );
    Eigen::Matrix< double, 3, 1> lagrangeCoefficients_acc_e2;
    lagrangeCoefficients_acc_e2 = lagrangeInterpolatorPointer_acc_e2->computeCoefficientsOfInterpolatingPolynomial();

    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_acc_e3 =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated_acc_e3 );
    Eigen::Matrix< double, 3, 1> lagrangeCoefficients_acc_e3;
    lagrangeCoefficients_acc_e3 = lagrangeInterpolatorPointer_acc_e3->computeCoefficientsOfInterpolatingPolynomial();

    // Compute powers of t to caculate current acceleration using the Lagrange coefficients
    double currentAcc_e1, currentAcc_e2, currentAcc_e3;
    // Matrix of current accelerations
    Eigen::Matrix< double, 3, 1 > t_powers;
    t_powers << currentTime*currentTime,
                currentTime,
                1.0;
    currentAcc_e1 = lagrangeCoefficients_acc_e1.adjoint() * t_powers;
    currentAcc_e2 = lagrangeCoefficients_acc_e2.adjoint() * t_powers;
    currentAcc_e3 = lagrangeCoefficients_acc_e3.adjoint() * t_powers;

    Eigen::Matrix< double, 1, 3 > currentAcc;
    currentAcc << currentAcc_e1, currentAcc_e2, currentAcc_e3;
    // Save current acceleration in spacecraft class for use in other functions
    spacecraftPointer_->setCurrentAcceleration( currentAcc );

    // Calculate derivatives
    Eigen::MatrixXd currentAccDerivative( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentAccDerivative.fill( 0.0 );
    // First derivative
    currentAccDerivative( 0, 0 ) = lagrangeCoefficients_acc_e1.adjoint() * dt_powers_1; // e1
    currentAccDerivative( 0, 1 ) = lagrangeCoefficients_acc_e2.adjoint() * dt_powers_1; // e2
    currentAccDerivative( 0, 2 ) = lagrangeCoefficients_acc_e3.adjoint() * dt_powers_1; // e3
    // Second derivative
    currentAccDerivative( 1, 0 ) = lagrangeCoefficients_acc_e1.adjoint() * dt_powers_2; // e1
    currentAccDerivative( 1, 1 ) = lagrangeCoefficients_acc_e2.adjoint() * dt_powers_2; // e2
    currentAccDerivative( 1, 2 ) = lagrangeCoefficients_acc_e3.adjoint() * dt_powers_2; // e3
    // Third to K-th derivatives are zero

    // Save current acceleration derivative in spacecraft for use in other functions
    spacecraftPointer_->setCurrentAccelerationDerivative( currentAccDerivative );
}

void DiscreteForceModelL3::updateCurrentForcesAndAccelerationsForTSI(
        double currentTime)
{
    /// Extract variables from classes
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix;
    thrustForceMatrix = spacecraftPointer_->getThrustForceMatrix();

    /// Specific relations for the specific low-thrust representation method
    // Find index in thrust profile that corresponds with the current time
    int index = tudat::basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                thrustForceMatrix.col( 0 ), currentTime);

    //??
//    std::cout << "index = " << index << std::endl;

    // Declarations
    double Time_prev; // Second-to-nearest left neighbor
    double Time_0; // Nearest left neighbor
    double Time_next; // Nearest right neighbor
    double Time_next2; // Second-to-nearest right neighbor. Added because 3rd order Lagrange requires 4 points.

    // First calculate previous, 0, next and next2 times in thrust profile to
    // simplify the next series of calculations
    Time_0 = thrustForceMatrix( index, 0 );
    // When index = 0 -> exception, use thrust value of first point twice
    if ( index == 0 )
    {
        // The difference between the first and second points is subtracted from first point and the result used as the new previous point
        Time_prev = Time_0 - ( thrustForceMatrix( index + 1, 0 ) - Time_0 );
    }
    else
    {
        Time_prev = thrustForceMatrix( index - 1, 0 ); // scalar of resulting thrust force
    }
    // Time at nearest right neighbour can always be calculated
    Time_next = thrustForceMatrix( index + 1, 0 );
    // When index = length( T ) -> exception, use last point twice
    if ( index > ( thrustForceMatrix.rows() - 2.5 ) )
    {
        // The difference between the last and second-to-last point is added to the last point and the result is used as the new next point
        Time_next2 = Time_next + ( Time_next - Time_0 );
    }
    else
    {
        Time_next2 = thrustForceMatrix( index + 2, 0 );
    }




    /////// Calculation of derivatives of resulting thrust at current time

    // Calculate derivative of the resulting current thrust force with the use
    // coefficients of the Lagrange polynomial

    // Declarations
    double T_res_prev;
    double T_res_0;
    double T_res_next;
    double T_res_next2; // Second-to-nearest right neighbor. Added because 3rd order Lagrange requires 4 points.

    T_res_0 = ( ( thrustForceMatrix.rightCols( 3 ) ).row( index ) ).norm();
    // When index == 0 -> exception, use first point twice
    if ( index == 0 )
    {
        T_res_prev = T_res_0; // The value of the previous point equals the value of the neares left neighbor (assumption)
    }
    else
    {
        T_res_prev = thrustForceMatrix.rightCols( 3 ).row( index - 1 ).norm(); // scalar of resulting thrust force
    }
    // Time at nearest right neighbour can always be calculated
    T_res_next = thrustForceMatrix.rightCols( 3 ).row( index + 1 ).norm();
    // When index = length( T ) -> exception, use last point twice
    if ( index > ( thrustForceMatrix.rows() - 2.5 )  )
    {
        T_res_next2 = T_res_next; // The last point is used twice ( assumption )
    }
    else
    {
        T_res_next2 = thrustForceMatrix.rightCols( 3 ).row( index + 2 ).norm();
    }

    // Calculate the coefficient of the Lagrange polynomial that goes through the
    // points of resulting thrust force
    Eigen::Matrix< double, 4, 2 > PointsToBeInterpolated;
    PointsToBeInterpolated << Time_prev,  T_res_prev,
                              Time_0,     T_res_0,
                              Time_next,  T_res_next,
                              Time_next2, T_res_next2;

    // Construct Lagrange Interpolator
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_T_res =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated );
    Eigen::Matrix< double, 4, 1> lagrangeCoefficients_T_res;
    lagrangeCoefficients_T_res = lagrangeInterpolatorPointer_T_res->computeCoefficientsOfInterpolatingPolynomial();

    // Calculate k-th derivative of the resulting thrust force at the current
    // time. The outcome will be used further on to calculate U8.
    int order = integrationSettingsPointer_->getOrderOfTaylorSeries();
    Eigen::MatrixXd dT_res( std::max( order, 3 ), 1 ); // At least 3 derivatives are calculated
    dT_res.fill( 0.0 );
    // Declarations
    Eigen::Matrix< double, 4, 1> t_powers;
    Eigen::Matrix< double, 4, 1> dt_powers_1;
    Eigen::Matrix< double, 4, 1> dt_powers_2;
    Eigen::Matrix< double, 4, 1> dt_powers_3;
    // Zeroth derivative = actual interpolated value of resulting thrust force
    t_powers <<    currentTime * currentTime * currentTime,
                   currentTime * currentTime,
                   currentTime,
                   1.0;
    double T_res_curr = lagrangeCoefficients_T_res.adjoint() * t_powers;
    // First derivative
    dt_powers_1 << 3.0 * currentTime * currentTime,
                   2.0 * currentTime,
                   1.0,
                   0.0;
    dT_res( 0, 0 ) = lagrangeCoefficients_T_res.adjoint() * dt_powers_1; // First derivative of resulting thrust force
    // Second derivative
    dt_powers_2 << 6.0 * currentTime,
                   2.0,
                   0.0,
                   0.0;
    dT_res( 1, 0 ) = lagrangeCoefficients_T_res.adjoint() * dt_powers_2; // Second derivative of resulting thrust force
    // Third derivative
    dt_powers_3 << 6.0,
                   0.0,
                   0.0,
                   0.0;
    dT_res( 2, 0 ) = lagrangeCoefficients_T_res.adjoint() * dt_powers_3; // Third derivative of resulting thrust force
    // Fourth to K-th derivatives are zero

    //??
//    std::cout << "t_powers = \n" << t_powers << std::endl;
//    std::cout << "dt_powers_1 = \n" << dt_powers_1 << std::endl;
//    std::cout << "dt_powers_2 = \n" << dt_powers_2 << std::endl;
//    std::cout << "dt_powers_3 = \n" << dt_powers_3 << std::endl;
//    std::cout << "dT_res = \n" << dT_res << std::endl;

    // Store derivatives of resultant low-thrust force in spacecraft class
    spacecraftPointer_->setCurrentResultantThrustForceDerivatives( dT_res );


    /// Calculation of derivatives of acceleration at current time

    // _curr means 'at current step of the integrator'
    // _0 means 'at nearest step in the thrust profile'

    // Extract mass from current state
    Eigen::MatrixXd currentState = spacecraftPointer_->getCurrentState();
    double currentMass = currentState( 7 );

    // Calculate mass flow at nearest left node in the thrust profile.
    // Assume that massDerivative_0 = currentMassDerivative (assumption) -> Assumption no longer needed!
    double currentMassDerivative;
    currentMassDerivative = - T_res_curr / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );

    // Declarations
    double mass_prev, mass_0, mass_next, mass_next2;
    // Calculate mass at previous and next steps of thrust profile
    mass_prev = currentMass + currentMassDerivative * ( Time_prev - currentTime );
    mass_0 = currentMass + currentMassDerivative * ( Time_0 - currentTime );
    mass_next = currentMass + currentMassDerivative * ( Time_next - currentTime );
    mass_next2 = currentMass + currentMassDerivative * ( Time_next2 - currentTime );

    //??
//    std::cout << "mass_prev = " << mass_prev << std::endl;
//    std::cout << "mass_0 = " << mass_0 << std::endl;
//    std::cout << "mass_next = " << mass_next << std::endl;
//    std::cout << "mass_next2 = " << mass_next2 << std::endl;

    // Transform thrust in velocity frame T_vf to thrust in USM frame
    double velocity_e1 = currentState( 8 );
    double velocity_e2 = currentState( 9 );

    // Calculate flight path angle gamma
    double sin_gamma, cos_gamma, gamma;
    sin_gamma = velocity_e1 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    cos_gamma = velocity_e2 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    gamma = std::atan2( sin_gamma, cos_gamma );

    // Calculate thrust vector at four consecutive nodes
    Eigen::Matrix< double, 1, 3 > thrust_USM_prev, thrust_USM_0, thrust_USM_next, thrust_USM_next2;
    // Create frameTransformation object
    FrameTransformation frameTransformation;
    thrust_USM_0 = frameTransformation.velocityFrameToUSMFrame( thrustForceMatrix.rightCols( 3 ).row( index ), gamma );
    // When index == 0 -> exception, use thrust value of first node twice (assumption)
    if ( index == 0 )
    {
        thrust_USM_prev = thrust_USM_0; // The value of the first node is used twice (assumption)
    }
    else
    {
        thrust_USM_prev = frameTransformation.velocityFrameToUSMFrame( thrustForceMatrix.rightCols( 3 ).row( index - 1 ), gamma );
    }
    // Thrust at nearest right neighbour can always be calculated
    thrust_USM_next = frameTransformation.velocityFrameToUSMFrame( thrustForceMatrix.rightCols( 3 ).row( index + 1 ), gamma );
    // When index == length( T ) -> exception, use thrust value of last node twice (assumption)
    if ( index > ( thrustForceMatrix.rows() - 2.5 ) )
    {
        thrust_USM_next2 = thrust_USM_next; // The thrust value of the last node is used twice (assumption)
    }
    else
    {
        thrust_USM_next2 = frameTransformation.velocityFrameToUSMFrame( thrustForceMatrix.rightCols( 3 ).row( index + 2 ), gamma );
    }

    // Calculation of the acceleration vectors at the previous, current and next steps of the thrust profile
    Eigen::Matrix< double, 3, 1 > acc_prev, acc_0, acc_next, acc_next2;
    acc_0 = thrust_USM_0 / mass_0; // 3x1 vector
    acc_prev = thrust_USM_prev / mass_prev; // 3x1 vector
    acc_next = thrust_USM_next / mass_next; // 3x1 vector
    acc_next2 = thrust_USM_next2 / mass_next2; // 3x1 vector

    //??
//    std::cout << "acc_prev = \n" << acc_prev << std::endl;
//    std::cout << "acc_0 = \n" << acc_0 << std::endl;
//    std::cout << "acc_next = \n" << acc_next << std::endl;
//    std::cout << "acc_next2 = \n" << acc_next2 << std::endl;

    // Calculate the coefficients of the Lagrange polynomial that goes
    // through the three points of the acceleration vector for each
    // direction (e1, e2 and e3).
    Eigen::Matrix< double, 4, 2 > PointsToBeInterpolated_acc_e1,
                                  PointsToBeInterpolated_acc_e2,
                                  PointsToBeInterpolated_acc_e3;
    PointsToBeInterpolated_acc_e1 << Time_prev,  acc_prev( 0),
                                     Time_0,     acc_0( 0 ),
                                     Time_next,  acc_next( 0 ),
                                     Time_next2, acc_next2( 0 );
    PointsToBeInterpolated_acc_e2 << Time_prev,  acc_prev( 1 ),
                                     Time_0,     acc_0( 1 ),
                                     Time_next,  acc_next( 1 ),
                                     Time_next2, acc_next2( 1 );
    PointsToBeInterpolated_acc_e3 << Time_prev,  acc_prev( 2 ),
                                     Time_0,     acc_0( 2 ),
                                     Time_next,  acc_next( 2 ),
                                     Time_next2, acc_next2( 2 );

    // Lagrange interpolation to calculate coefficients
    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_acc_e1 =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated_acc_e1 );
    Eigen::Matrix< double, 4, 1> lagrangeCoefficients_acc_e1;
    lagrangeCoefficients_acc_e1 = lagrangeInterpolatorPointer_acc_e1->computeCoefficientsOfInterpolatingPolynomial();

    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_acc_e2 =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated_acc_e2 );
    Eigen::Matrix< double, 4, 1> lagrangeCoefficients_acc_e2;
    lagrangeCoefficients_acc_e2 = lagrangeInterpolatorPointer_acc_e2->computeCoefficientsOfInterpolatingPolynomial();

    LagrangeInterpolatorPointer lagrangeInterpolatorPointer_acc_e3 =
            boost::make_shared< LagrangeInterpolator >( PointsToBeInterpolated_acc_e3 );
    Eigen::Matrix< double, 4, 1> lagrangeCoefficients_acc_e3;
    lagrangeCoefficients_acc_e3 = lagrangeInterpolatorPointer_acc_e3->computeCoefficientsOfInterpolatingPolynomial();

    // Compute powers of t to caculate current acceleration using the Lagrange coefficients
    double currentAcc_e1, currentAcc_e2, currentAcc_e3;
    // Matrix of current accelerations
    currentAcc_e1 = lagrangeCoefficients_acc_e1.adjoint() * t_powers;
    currentAcc_e2 = lagrangeCoefficients_acc_e2.adjoint() * t_powers;
    currentAcc_e3 = lagrangeCoefficients_acc_e3.adjoint() * t_powers;

    Eigen::Matrix< double, 1, 3 > currentAcc;
    currentAcc << currentAcc_e1, currentAcc_e2, currentAcc_e3;

    //??
//    std::cout << "PointsToBeInterpolated_acc_e1 = \n" << PointsToBeInterpolated_acc_e1 << std::endl;
//    std::cout << "PointsToBeInterpolated_acc_e2 = \n" << PointsToBeInterpolated_acc_e2 << std::endl;
//    std::cout << "PointsToBeInterpolated_acc_e3 = \n" << PointsToBeInterpolated_acc_e3 << std::endl;
//    std::cout << "lagrangeCoefficients_acc_e1.adjoint = \n" << lagrangeCoefficients_acc_e1.adjoint() << std::endl;
//    std::cout << "lagrangeCoefficients_acc_e2.adjoint = \n" << lagrangeCoefficients_acc_e2.adjoint() << std::endl;
//    std::cout << "lagrangeCoefficients_acc_e3.adjoint = \n" << lagrangeCoefficients_acc_e3.adjoint() << std::endl;
//    std::cout << "currentAcc = \n" << currentAcc << std::endl;

    // Save current acceleration in spacecraft class for use in other functions
    spacecraftPointer_->setCurrentAcceleration( currentAcc );

    // Calculate derivatives
    Eigen::MatrixXd currentAccDerivative( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentAccDerivative.fill( 0.0 );
    // First derivative
    currentAccDerivative( 0, 0 ) = lagrangeCoefficients_acc_e1.adjoint() * dt_powers_1; // e1
    currentAccDerivative( 0, 1 ) = lagrangeCoefficients_acc_e2.adjoint() * dt_powers_1; // e2
    currentAccDerivative( 0, 2 ) = lagrangeCoefficients_acc_e3.adjoint() * dt_powers_1; // e3
    // Second derivative
    currentAccDerivative( 1, 0 ) = lagrangeCoefficients_acc_e1.adjoint() * dt_powers_2; // e1
    currentAccDerivative( 1, 1 ) = lagrangeCoefficients_acc_e2.adjoint() * dt_powers_2; // e2
    currentAccDerivative( 1, 2 ) = lagrangeCoefficients_acc_e3.adjoint() * dt_powers_2; // e3
    // Third derivative
    currentAccDerivative( 2, 0 ) = lagrangeCoefficients_acc_e1.adjoint() * dt_powers_3; // e1
    currentAccDerivative( 2, 1 ) = lagrangeCoefficients_acc_e2.adjoint() * dt_powers_3; // e2
    currentAccDerivative( 2, 2 ) = lagrangeCoefficients_acc_e3.adjoint() * dt_powers_3; // e3
    // Fourth to K-th derivatives are zero

    // Save current acceleration derivative in spacecraft for use in other functions
    spacecraftPointer_->setCurrentAccelerationDerivative( currentAccDerivative );
}

void DiscreteForceModelCSI::updateCurrentForcesAndAccelerationsForTSI(
        double currentTime)
{

    //??
    // std::cout << "label 1" << std::endl;

    /// Extract variables from classes
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix;
    thrustForceMatrix = spacecraftPointer_->getThrustForceMatrix();

    /// Specific relations for the specific low-thrust representation method
    // Find index in thrust profile that corresponds with the current time
    int index = tudat::basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                thrustForceMatrix.col( 0 ), currentTime);
    double Time_0 = thrustForceMatrix( index, 0 ); // calculate time of left node

    //??
    // std::cout << "label 2" << std::endl;

    // Get CSI coefficients of resultant and x, y and z components in velocity frame from spacecraft pointer
    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX,
                                         CSIcoefficientsOfThrustForceInVelocityFrameY,
                                         CSIcoefficientsOfThrustForceInVelocityFrameZ,
                                         CSIcoefficientsOfResultantThrustForce;

    CSIcoefficientsOfThrustForceInVelocityFrameX =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameX(  );
    CSIcoefficientsOfThrustForceInVelocityFrameY =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameY(  );
    CSIcoefficientsOfThrustForceInVelocityFrameZ =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameZ(  );
    CSIcoefficientsOfResultantThrustForce =
            spacecraftPointer_->getCSIcoefficientsOfResultantThrustForce(  );

    /// Calculate values of resultant and x, y and z components in velocity frame
    double time_temp = currentTime - Time_0;
    double currentThrustForceInVelocityFrameX =
              CSIcoefficientsOfThrustForceInVelocityFrameX[ 0 ][ index ]
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 1 ][ index ] * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ] * time_temp * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameY =
            CSIcoefficientsOfThrustForceInVelocityFrameY[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameZ =
            CSIcoefficientsOfThrustForceInVelocityFrameZ[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp * time_temp * time_temp;
//    double currentResultantThrustForce =
//            CSIcoefficientsOfResultantThrustForce[ 0 ][ index ]
//          + CSIcoefficientsOfResultantThrustForce[ 1 ][ index ] * time_temp
//          + CSIcoefficientsOfResultantThrustForce[ 2 ][ index ] * time_temp * time_temp
//          + CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentResultantThrustForce = std::sqrt(
                currentThrustForceInVelocityFrameX * currentThrustForceInVelocityFrameX +
                currentThrustForceInVelocityFrameY * currentThrustForceInVelocityFrameY +
                currentThrustForceInVelocityFrameZ * currentThrustForceInVelocityFrameZ );

    //??
//    if ( currentResultantThrustForce < 0.0 )
//    {
//        std::cout<< currentResultantThrustForce << std::endl;
//    }

    // Save current resultant thrust force in spacecraft pointer
    spacecraftPointer_->setCurrentResultantThrustForce( currentResultantThrustForce );

    /// Calculate the current derivatives of the resultant thrust force
    int order = integrationSettingsPointer_->getOrderOfTaylorSeries();
    Eigen::MatrixXd currentDerivativesOfResultantThrustForce( std::max( order, 3 ), 1 ); // At least 3 derivatives are calculated
    currentDerivativesOfResultantThrustForce.fill( 0.0 );

    currentDerivativesOfResultantThrustForce( 0, 0 ) =
                    CSIcoefficientsOfResultantThrustForce[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfResultantThrustForce[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfResultantThrustForce( 1, 0 ) =
              2.0 * CSIcoefficientsOfResultantThrustForce[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp;
    currentDerivativesOfResultantThrustForce( 2, 0 ) =
              6.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ];

    // Store derivatives of resultant low-thrust force in spacecraft class
    spacecraftPointer_->setCurrentResultantThrustForceDerivatives( currentDerivativesOfResultantThrustForce );

    //??
    // std::cout << "label 3" << std::endl;


    /// Calculation of the current thrust acceleration in USM frame components

    // Extract mass from current state
    Eigen::MatrixXd currentState = spacecraftPointer_->getCurrentState();
    double currentMass = currentState( 7 );

    // Calculate flight path angle gamma required for transformation from VF to USM frame
    double velocity_e1 = currentState( 8 );
    double velocity_e2 = currentState( 9 );
    double sin_gamma, cos_gamma, gamma;
    sin_gamma = velocity_e1 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    cos_gamma = velocity_e2 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    gamma = std::atan2( sin_gamma, cos_gamma );

    // Create frameTransformation object
    FrameTransformation frameTransformation;

    // Create vector of currentThrustForceInVelocityFrame
    Eigen::Matrix< double, 1, 3 > currentThrustForceInVelocityFrame;
    currentThrustForceInVelocityFrame << currentThrustForceInVelocityFrameX,
                                         currentThrustForceInVelocityFrameY,
                                         currentThrustForceInVelocityFrameZ;

    // Compute current thrust force in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustForceInUSMframe =
            frameTransformation.velocityFrameToUSMFrame( currentThrustForceInVelocityFrame, gamma );

    // Compute current thrust acceleration in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustAccelerationInUSMframe =
            currentThrustForceInUSMframe / currentMass;

    // Save current thrust acceleration in USM frame components in spacecraft pointer
    spacecraftPointer_->setCurrentAcceleration( currentThrustAccelerationInUSMframe );



    //??
    // std::cout << "label 4" << std::endl;

    /// Calculation of current derivatives of the thrust accelerations in USM frame components

    // Calculate the current derivatives of the thrust force components in the velocity frame
    Eigen::MatrixXd currentDerivativesOfThrustInVelocityFrame( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustInVelocityFrame.fill( 0.0 );
    // x-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 0 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameX[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 0 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 0 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ];
    // y-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 1 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameY[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 1 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 1 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ];
    // z-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 2 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameZ[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 2 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 2 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ];

    //??
    // std::cout << "label 5" << std::endl;

    // Transform the current derivatives of the thrust force components to USM frame
    Eigen::MatrixXd currentDerivativesOfThrustInUSMframe( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustInUSMframe.fill( 0.0 );
    Eigen::Matrix< double, 1, 3 > temp_matrix;
    for ( int i = 0; i < 3; i++ ) // because there are only 3 non-zero derivatives
    {
        // Assumption that gamma_dot = 0 for near circular orbits of GTOC3
        temp_matrix = frameTransformation.velocityFrameToUSMFrame(
                    currentDerivativesOfThrustInVelocityFrame.row( i ), gamma );
        currentDerivativesOfThrustInUSMframe.row( i ) = temp_matrix;
    }

    //??
    // std::cout << "label 6" << std::endl;

    // Calculate derivatives of current thrust acceleration in USM frame components
    Eigen::VectorXd currentMassDerivatives( 3 ); // first three derivatives of mass
    currentMassDerivatives( 0 ) = - currentResultantThrustForce
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );
    currentMassDerivatives( 1 ) = - currentDerivativesOfResultantThrustForce( 0, 0 )
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );
    currentMassDerivatives( 2 ) = - currentDerivativesOfResultantThrustForce( 1, 0 )
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );

    //??
    // // std::cout << "label 7" << std::endl;

    Eigen::MatrixXd currentDerivativesOfThrustAccelerationInUSMframe( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustAccelerationInUSMframe.fill( 0.0 );

    //??
    // // std::cout << "label 8" << std::endl;

    for ( int j = 0; j < 3; j++ )
    {
        // First derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 0, j ) =
                ( currentDerivativesOfThrustInUSMframe( 0, j ) * currentMass
                  - currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass );
        // Second derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 1, j ) =
                ( currentDerivativesOfThrustInUSMframe( 1, j ) * currentMass * currentMass
                  - currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMass
                  - 2.0 * currentMass * currentMassDerivatives( 0 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                  + 2.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass * currentMass );
        // Third derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 2, j ) =
                ( ( currentDerivativesOfThrustInUSMframe( 2, j ) * currentMass * currentMass
                    - currentMassDerivatives( 2 ) * currentThrustForceInUSMframe( j ) * currentMass
                    - currentMassDerivatives( 1 ) * currentDerivativesOfThrustInUSMframe( 0, j ) * currentMass
                    - currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMassDerivatives( 0 )
                    - 2.0 * currentMass * currentMassDerivatives( 1 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                    + 4.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) )
                  * currentMass
                  - 3.0 * currentDerivativesOfThrustInUSMframe( 1, j ) * currentMass * currentMass
                  + 3.0 * currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMass
                  + 6.0 * currentMass * currentMassDerivatives( 0 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                  - 6.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass * currentMass * currentMass );
    }


    // Save current derivatives of thrust acceleration in USM frame components in spacecraft pointer
    spacecraftPointer_->setCurrentAccelerationDerivative( currentDerivativesOfThrustAccelerationInUSMframe );

}

void DiscreteForceModelCSIcorrectGamma::updateCurrentForcesAndAccelerationsForTSI(
        double currentTime)
{

    //??
    // std::cout << "label 1" << std::endl;

    /// Extract variables from classes
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix;
    thrustForceMatrix = spacecraftPointer_->getThrustForceMatrix();

    /// Specific relations for the specific low-thrust representation method
    // Find index in thrust profile that corresponds with the current time
    int index = tudat::basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                thrustForceMatrix.col( 0 ), currentTime);
    double Time_0 = thrustForceMatrix( index, 0 ); // calculate time of left node

    //??
    // std::cout << "label 2" << std::endl;

    // Get CSI coefficients of resultant and x, y and z components in velocity frame from spacecraft pointer
    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX,
                                         CSIcoefficientsOfThrustForceInVelocityFrameY,
                                         CSIcoefficientsOfThrustForceInVelocityFrameZ,
                                         CSIcoefficientsOfResultantThrustForce;

    CSIcoefficientsOfThrustForceInVelocityFrameX =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameX(  );
    CSIcoefficientsOfThrustForceInVelocityFrameY =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameY(  );
    CSIcoefficientsOfThrustForceInVelocityFrameZ =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameZ(  );
    CSIcoefficientsOfResultantThrustForce =
            spacecraftPointer_->getCSIcoefficientsOfResultantThrustForce(  );

    /// Calculate values of resultant and x, y and z components in velocity frame
    double time_temp = currentTime - Time_0;
    double currentThrustForceInVelocityFrameX =
              CSIcoefficientsOfThrustForceInVelocityFrameX[ 0 ][ index ]
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 1 ][ index ] * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ] * time_temp * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameY =
            CSIcoefficientsOfThrustForceInVelocityFrameY[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameZ =
            CSIcoefficientsOfThrustForceInVelocityFrameZ[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp * time_temp * time_temp;
    //    double currentResultantThrustForce =
    //            CSIcoefficientsOfResultantThrustForce[ 0 ][ index ]
    //          + CSIcoefficientsOfResultantThrustForce[ 1 ][ index ] * time_temp
    //          + CSIcoefficientsOfResultantThrustForce[ 2 ][ index ] * time_temp * time_temp
    //          + CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp * time_temp * time_temp;
        double currentResultantThrustForce = std::sqrt(
                    currentThrustForceInVelocityFrameX * currentThrustForceInVelocityFrameX +
                    currentThrustForceInVelocityFrameY * currentThrustForceInVelocityFrameY +
                    currentThrustForceInVelocityFrameZ * currentThrustForceInVelocityFrameZ );

    //??
    if ( currentResultantThrustForce < 0.0 )
    {
        std::cout<< currentResultantThrustForce << std::endl;
    }

    // Save current resultant thrust force in spacecraft pointer
    spacecraftPointer_->setCurrentResultantThrustForce( currentResultantThrustForce );

    /// Calculate the current derivatives of the resultant thrust force
    int order = integrationSettingsPointer_->getOrderOfTaylorSeries();
    Eigen::MatrixXd currentDerivativesOfResultantThrustForce( std::max( order, 3 ), 1 ); // At least 3 derivatives are calculated
    currentDerivativesOfResultantThrustForce.fill( 0.0 );

    currentDerivativesOfResultantThrustForce( 0, 0 ) =
                    CSIcoefficientsOfResultantThrustForce[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfResultantThrustForce[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfResultantThrustForce( 1, 0 ) =
              2.0 * CSIcoefficientsOfResultantThrustForce[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp;
    currentDerivativesOfResultantThrustForce( 2, 0 ) =
              6.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ];

    // Store derivatives of resultant low-thrust force in spacecraft class
    spacecraftPointer_->setCurrentResultantThrustForceDerivatives( currentDerivativesOfResultantThrustForce );

    //??
    // std::cout << "label 3" << std::endl;


    /// Calculation of the current thrust acceleration in USM frame components

    // Extract mass from current state
    Eigen::MatrixXd currentState = spacecraftPointer_->getCurrentState();
    double currentMass = currentState( 7 );

    // Calculate flight path angle gamma required for transformation from VF to USM frame
    double velocity_e1 = currentState( 8 );
    double velocity_e2 = currentState( 9 );
    double sin_gamma, cos_gamma, gamma;
    sin_gamma = velocity_e1 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    cos_gamma = velocity_e2 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    gamma = std::atan2( sin_gamma, cos_gamma );

    // Create frameTransformation object
    FrameTransformation frameTransformation;

    // Create vector of currentThrustForceInVelocityFrame
    Eigen::Matrix< double, 1, 3 > currentThrustForceInVelocityFrame;
    currentThrustForceInVelocityFrame << currentThrustForceInVelocityFrameX,
                                         currentThrustForceInVelocityFrameY,
                                         currentThrustForceInVelocityFrameZ;

    // Compute current thrust force in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustForceInUSMframe =
            frameTransformation.velocityFrameToUSMFrame( currentThrustForceInVelocityFrame, gamma );

    // Compute current thrust acceleration in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustAccelerationInUSMframe =
            currentThrustForceInUSMframe / currentMass;

    // Save current thrust acceleration in USM frame components in spacecraft pointer
    spacecraftPointer_->setCurrentAcceleration( currentThrustAccelerationInUSMframe );



    //??
    // std::cout << "label 4" << std::endl;

    /// Calculation of current derivatives of the thrust accelerations in USM frame components

    // Calculate the current derivatives of the thrust force components in the velocity frame
    Eigen::MatrixXd currentDerivativesOfThrustInVelocityFrame( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustInVelocityFrame.fill( 0.0 );
    // x-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 0 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameX[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 0 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 0 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ];
    // y-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 1 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameY[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 1 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 1 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ];
    // z-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 2 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameZ[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 2 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 2 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ];

    //??
    // std::cout << "label 5" << std::endl;

    // Transform the current derivatives of the thrust force components to USM frame
    Eigen::MatrixXd currentDerivativesOfThrustInUSMframe( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustInUSMframe.fill( 0.0 );
    Eigen::Matrix< double, 1, 3 > temp_matrix;

    // Extract derivative matrix of current coefficient matrix
    Eigen::MatrixXd currentCoefficientMatrix = spacecraftPointer_->getCurrentCoefficientMatrix();
    double v_r = velocity_e1;
    double dv_r = currentCoefficientMatrix( 1, 8 ); // 0th row is current state, 1st row is first derivative of state
    double ddv_r = currentCoefficientMatrix( 2, 8 );
    double dddv_r = currentCoefficientMatrix( 3, 8 );
    double v_tr = velocity_e2;
    double dv_tr = currentCoefficientMatrix( 1, 9 ); // 0th row is current state, 1st row is first derivative of state
    double ddv_tr = currentCoefficientMatrix( 2, 9 );
    double dddv_tr = currentCoefficientMatrix( 3, 9 );

    // Analytic calculation of the first three derivatives of the flight path angle gamma
    double gammaDot = ( dv_r /v_tr - (v_r*dv_tr)/(v_tr * v_tr))/((v_r * v_r)/(v_tr * v_tr) + 1.0 );
    double gammaDDot = (ddv_r/v_tr - (2.0*dv_r*dv_tr)/(v_tr * v_tr)
                        + (2.0*v_r*dv_tr * dv_tr)/(v_tr * v_tr * v_tr) - (v_r*ddv_tr)/(v_tr * v_tr))/((v_r * v_r)/(v_tr * v_tr) + 1.0)
            - ((dv_r/v_tr - (v_r*dv_tr)/(v_tr * v_tr))
               *((2.0*v_r*dv_r)/(v_tr * v_tr) - (2.0*(v_r * v_r)*dv_tr)/(v_tr * v_tr * v_tr)))/
            std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 2.0 );
    double gammaDDDot = (2.0*(dv_r/v_tr
                            - (v_r*dv_tr)/(v_tr * v_tr))*std::pow( ((2.0*v_r*dv_r)/(v_tr * v_tr)
                                                                    - (2.0*(v_r * v_r)*dv_tr)/(v_tr * v_tr * v_tr)), 2.0) )/
            std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 3.0)
            - (2.0*((2.0*v_r*dv_r)/(v_tr * v_tr) - (2.0*(v_r * v_r)*dv_tr)/(v_tr * v_tr * v_tr))*(ddv_r/v_tr
                                                                               - (2.0*dv_r*dv_tr)/(v_tr * v_tr)
                                                                               + (2.0*v_r*dv_tr * dv_tr)/(v_tr * v_tr * v_tr)
                                                                               - (v_r*ddv_tr)/(v_tr * v_tr)))/
            std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 2.0 )
            - ((dv_r/v_tr - (v_r*dv_tr)/(v_tr * v_tr))*((2.0*dv_r * dv_r)/(v_tr * v_tr)
                                                    + (6.0*(v_r * v_r)*dv_tr * dv_tr)/(v_tr * v_tr * v_tr * v_tr)
                                                    - (2.0*(v_r * v_r)*ddv_tr)/(v_tr * v_tr * v_tr)
                                                    + (2.0*v_r*ddv_r)/(v_tr * v_tr)
                                                    - (8.0*v_r*dv_r*dv_tr)/(v_tr * v_tr * v_tr)))
            /std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 2.0 ) - ((v_r*dddv_tr)/(v_tr * v_tr) - dddv_r/v_tr
                                           - (6.0*dv_r*dv_tr * dv_tr)/(v_tr * v_tr * v_tr) + (3.0*dv_r*ddv_tr)/(v_tr * v_tr)
                                           + (3.0*dv_tr*ddv_r)/(v_tr * v_tr) + (6.0*v_r*dv_tr * dv_tr * dv_tr)/(v_tr * v_tr * v_tr * v_tr)
                                           - (6.0*v_r*dv_tr*ddv_tr)/(v_tr * v_tr * v_tr))/((v_r * v_r)/(v_tr * v_tr) + 1.0);

    Eigen::Matrix< double, 3, 3> transformationMatrixVftoUSMframe;
    transformationMatrixVftoUSMframe <<
            std::sin( gamma ), -std::cos( gamma ), 0.0,
            std::cos( gamma ), std::sin( gamma ), 0.0,
            0.0, 0.0, 1.0;
    Eigen::Matrix< double, 3, 3> firstDerivativeOfTransformationMatrixVfToUSMframe;
    firstDerivativeOfTransformationMatrixVfToUSMframe <<
            gammaDot * std::cos( gamma ), gammaDot * std::sin( gamma ), 0.0,
            -gammaDot * sin( gamma ), gammaDot * std::cos( gamma ), 0.0,
            0.0, 0.0, 0.0;
    Eigen::Matrix< double, 3, 3> secondDerivativeOfTransformationMatrixVfToUSMframe;
    secondDerivativeOfTransformationMatrixVfToUSMframe <<
            gammaDDot * std::cos( gamma ) - gammaDot * gammaDot * std::sin( gamma ),
            gammaDDot * std::sin( gamma ) + gammaDot * gammaDot * std::cos( gamma ), 0.0,
            -gammaDDot * std::sin( gamma ) - gammaDot * gammaDot * std::cos( gamma ),
            gammaDDot * std::cos( gamma ) - gammaDot * gammaDot * std::sin( gamma ), 0.0,
            0.0, 0.0, 0.0;
    Eigen::Matrix< double, 3, 3> thirdDerivativeOfTransformationMatrixVfToUSMframe;
    thirdDerivativeOfTransformationMatrixVfToUSMframe <<
            gammaDDDot * std::cos( gamma )
            - 3.0 * gammaDDot * gammaDot * std::sin( gamma )
            - gammaDot * gammaDot * gammaDot * std::cos( gamma ),
            gammaDDDot * std::sin( gamma )
            + 3.0 * gammaDDot * gammaDot * std::cos( gamma )
            - gammaDot * gammaDot * gammaDot * std::sin( gamma ), 0.0,
            - gammaDDDot * std::sin( gamma )
            - 3.0 * gammaDDot * gammaDot * std::cos( gamma )
            + gammaDot * gammaDot * gammaDot * std::sin( gamma ),
            gammaDDDot * std::cos( gamma )
            - 3.0 * gammaDDot * gammaDot * std::sin( gamma )
            - gammaDot * gammaDot * gammaDot * std::cos( gamma ), 0.0,
            0.0, 0.0, 0.0;


    // First derivative
    currentDerivativesOfThrustInUSMframe.row( 0 ) =
            ( firstDerivativeOfTransformationMatrixVfToUSMframe * currentThrustForceInVelocityFrame.adjoint( )
              + transformationMatrixVftoUSMframe * currentDerivativesOfThrustInVelocityFrame.row( 0 ).adjoint( ) ).adjoint( );
    currentDerivativesOfThrustInUSMframe.row( 1 ) =
            ( secondDerivativeOfTransformationMatrixVfToUSMframe * currentThrustForceInVelocityFrame.adjoint( )
              + 2.0 * firstDerivativeOfTransformationMatrixVfToUSMframe *
                      currentDerivativesOfThrustInVelocityFrame.row( 0 ).adjoint( )
              + transformationMatrixVftoUSMframe * currentDerivativesOfThrustInVelocityFrame.row( 1 ).adjoint( ) ).adjoint( );
    currentDerivativesOfThrustInUSMframe.row( 2 ) =
            ( thirdDerivativeOfTransformationMatrixVfToUSMframe * currentThrustForceInVelocityFrame.adjoint( )
              + 3.0 * secondDerivativeOfTransformationMatrixVfToUSMframe *
                      currentDerivativesOfThrustInVelocityFrame.row( 0 ).adjoint( )
              + 3.0 * firstDerivativeOfTransformationMatrixVfToUSMframe *
                      currentDerivativesOfThrustInVelocityFrame.row( 1 ).adjoint( )
              + transformationMatrixVftoUSMframe * currentDerivativesOfThrustInVelocityFrame.row( 2 ).adjoint( ) ).adjoint( );

    //??
    // std::cout << "label 6" << std::endl;

    // Calculate derivatives of current thrust acceleration in USM frame components
    Eigen::VectorXd currentMassDerivatives( 3 ); // first three derivatives of mass
    currentMassDerivatives( 0 ) = - currentResultantThrustForce
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );
    currentMassDerivatives( 1 ) = - currentDerivativesOfResultantThrustForce( 0, 0 )
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );
    currentMassDerivatives( 2 ) = - currentDerivativesOfResultantThrustForce( 1, 0 )
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );

    //??
    // // std::cout << "label 7" << std::endl;

    Eigen::MatrixXd currentDerivativesOfThrustAccelerationInUSMframe( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustAccelerationInUSMframe.fill( 0.0 );

    //??
    // // std::cout << "label 8" << std::endl;

    for ( int j = 0; j < 3; j++ )
    {
        // First derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 0, j ) =
                ( currentDerivativesOfThrustInUSMframe( 0, j ) * currentMass
                  - currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass );
        // Second derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 1, j ) =
                ( currentDerivativesOfThrustInUSMframe( 1, j ) * currentMass * currentMass
                  - currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMass
                  - 2.0 * currentMass * currentMassDerivatives( 0 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                  + 2.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass * currentMass );
        // Third derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 2, j ) =
                ( ( currentDerivativesOfThrustInUSMframe( 2, j ) * currentMass * currentMass
                    - currentMassDerivatives( 2 ) * currentThrustForceInUSMframe( j ) * currentMass
                    - currentMassDerivatives( 1 ) * currentDerivativesOfThrustInUSMframe( 0, j ) * currentMass
                    - currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMassDerivatives( 0 )
                    - 2.0 * currentMass * currentMassDerivatives( 1 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                    + 4.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) )
                  * currentMass
                  - 3.0 * currentDerivativesOfThrustInUSMframe( 1, j ) * currentMass * currentMass
                  + 3.0 * currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMass
                  + 6.0 * currentMass * currentMassDerivatives( 0 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                  - 6.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass * currentMass * currentMass );
    }


    // Save current derivatives of thrust acceleration in USM frame components in spacecraft pointer
    spacecraftPointer_->setCurrentAccelerationDerivative( currentDerivativesOfThrustAccelerationInUSMframe );

}



void DiscreteForceModelCSIimprovedMass::updateCurrentForcesAndAccelerationsForTSI(
        double currentTime)
{

    //??
    // std::cout << "label 1" << std::endl;

    /// Extract variables from classes
    Eigen::Matrix< double, Eigen::Dynamic, 4 > thrustForceMatrix;
    thrustForceMatrix = spacecraftPointer_->getThrustForceMatrix();

    /// Specific relations for the specific low-thrust representation method
    // Find index in thrust profile that corresponds with the current time
    int index = tudat::basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                thrustForceMatrix.col( 0 ), currentTime);
    double Time_0 = thrustForceMatrix( index, 0 ); // calculate time of left node

    //??
    // std::cout << "label 2" << std::endl;

    // Get CSI coefficients of resultant and x, y and z components in velocity frame from spacecraft pointer
    std::vector< std::vector< double > > CSIcoefficientsOfThrustForceInVelocityFrameX,
                                         CSIcoefficientsOfThrustForceInVelocityFrameY,
                                         CSIcoefficientsOfThrustForceInVelocityFrameZ,
                                         CSIcoefficientsOfResultantThrustForce;

    CSIcoefficientsOfThrustForceInVelocityFrameX =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameX(  );
    CSIcoefficientsOfThrustForceInVelocityFrameY =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameY(  );
    CSIcoefficientsOfThrustForceInVelocityFrameZ =
            spacecraftPointer_->getCSIcoefficientsOfThrustForceInVelocityFrameZ(  );
    CSIcoefficientsOfResultantThrustForce =
            spacecraftPointer_->getCSIcoefficientsOfResultantThrustForce(  );

    /// Calculate values of resultant and x, y and z components in velocity frame
    double time_temp = currentTime - Time_0;
    double currentThrustForceInVelocityFrameX =
              CSIcoefficientsOfThrustForceInVelocityFrameX[ 0 ][ index ]
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 1 ][ index ] * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ] * time_temp * time_temp
            + CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameY =
            CSIcoefficientsOfThrustForceInVelocityFrameY[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp * time_temp * time_temp;
    double currentThrustForceInVelocityFrameZ =
            CSIcoefficientsOfThrustForceInVelocityFrameZ[ 0 ][ index ]
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 1 ][ index ] * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ] * time_temp * time_temp
          + CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp * time_temp * time_temp;
    //    double currentResultantThrustForce =
    //            CSIcoefficientsOfResultantThrustForce[ 0 ][ index ]
    //          + CSIcoefficientsOfResultantThrustForce[ 1 ][ index ] * time_temp
    //          + CSIcoefficientsOfResultantThrustForce[ 2 ][ index ] * time_temp * time_temp
    //          + CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp * time_temp * time_temp;
        double currentResultantThrustForce = std::sqrt(
                    currentThrustForceInVelocityFrameX * currentThrustForceInVelocityFrameX +
                    currentThrustForceInVelocityFrameY * currentThrustForceInVelocityFrameY +
                    currentThrustForceInVelocityFrameZ * currentThrustForceInVelocityFrameZ );

    //??
    if ( currentResultantThrustForce < 0.0 )
    {
        std::cout<< currentResultantThrustForce << std::endl;
    }

    // Save current resultant thrust force in spacecraft pointer
    spacecraftPointer_->setCurrentResultantThrustForce( currentResultantThrustForce );

    /// Calculate the current derivatives of the resultant thrust force
    int order = integrationSettingsPointer_->getOrderOfTaylorSeries();
    Eigen::MatrixXd currentDerivativesOfResultantThrustForce( std::max( order, 3 ), 1 ); // At least 3 derivatives are calculated
    currentDerivativesOfResultantThrustForce.fill( 0.0 );

    // Here the current derivatives are calculated using h (currentAlternativeTime) instead of t (currentTime) as the time variable,
    // with h = t / dt resulting in (h - h_0)^k = 1 for fixed step sizes. This ensures that the calculation of the Taylor series
    // for the mass is less prune to inaccuracies.
    double stepSize = integrationSettingsPointer_->getStepSize();

    currentDerivativesOfResultantThrustForce( 0, 0 ) =
            stepSize * (
                    CSIcoefficientsOfResultantThrustForce[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfResultantThrustForce[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp * time_temp
            );
    currentDerivativesOfResultantThrustForce( 1, 0 ) =
            stepSize * stepSize * (
              2.0 * CSIcoefficientsOfResultantThrustForce[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ] * time_temp
            );
    currentDerivativesOfResultantThrustForce( 2, 0 ) =
            stepSize * stepSize * stepSize *
              6.0 * CSIcoefficientsOfResultantThrustForce[ 3 ][ index ];

    // Store derivatives of resultant low-thrust force in spacecraft class
    spacecraftPointer_->setCurrentResultantThrustForceDerivatives( currentDerivativesOfResultantThrustForce );

    //??
    // std::cout << "label 3" << std::endl;


    /// Calculation of the current thrust acceleration in USM frame components

    // Extract mass from current state
    Eigen::MatrixXd currentState = spacecraftPointer_->getCurrentState();
    double currentMass = currentState( 7 );

    // Calculate flight path angle gamma required for transformation from VF to USM frame
    double velocity_e1 = currentState( 8 );
    double velocity_e2 = currentState( 9 );
    double sin_gamma, cos_gamma, gamma;
    sin_gamma = velocity_e1 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    cos_gamma = velocity_e2 / std::sqrt( velocity_e1 * velocity_e1 + velocity_e2 * velocity_e2 );
    gamma = std::atan2( sin_gamma, cos_gamma );

    // Create frameTransformation object
    FrameTransformation frameTransformation;

    // Create vector of currentThrustForceInVelocityFrame
    Eigen::Matrix< double, 1, 3 > currentThrustForceInVelocityFrame;
    currentThrustForceInVelocityFrame << currentThrustForceInVelocityFrameX,
                                         currentThrustForceInVelocityFrameY,
                                         currentThrustForceInVelocityFrameZ;

    // Compute current thrust force in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustForceInUSMframe =
            frameTransformation.velocityFrameToUSMFrame( currentThrustForceInVelocityFrame, gamma );

    // Compute current thrust acceleration in USM frame components
    Eigen::Matrix< double, 1, 3 > currentThrustAccelerationInUSMframe =
            currentThrustForceInUSMframe / currentMass;

    // Save current thrust acceleration in USM frame components in spacecraft pointer
    spacecraftPointer_->setCurrentAcceleration( currentThrustAccelerationInUSMframe );



    //??
    // std::cout << "label 4" << std::endl;

    /// Calculation of current derivatives of the thrust accelerations in USM frame components

    // Calculate the current derivatives of the thrust force components in the velocity frame
    Eigen::MatrixXd currentDerivativesOfThrustInVelocityFrame( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustInVelocityFrame.fill( 0.0 );
    // x-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 0 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameX[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 0 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 0 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameX[ 3 ][ index ];
    // y-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 1 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameY[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 1 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 1 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameY[ 3 ][ index ];
    // z-direction
    currentDerivativesOfThrustInVelocityFrame( 0, 2 ) =
                    CSIcoefficientsOfThrustForceInVelocityFrameZ[ 1 ][ index ]
            + 2.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ] * time_temp
            + 3.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 1, 2 ) =
              2.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 2 ][ index ]
            + 6.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ] * time_temp;
    currentDerivativesOfThrustInVelocityFrame( 2, 2 ) =
              6.0 * CSIcoefficientsOfThrustForceInVelocityFrameZ[ 3 ][ index ];

    //??
    // std::cout << "label 5" << std::endl;

    // Transform the current derivatives of the thrust force components to USM frame
    Eigen::MatrixXd currentDerivativesOfThrustInUSMframe( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustInUSMframe.fill( 0.0 );
    Eigen::Matrix< double, 1, 3 > temp_matrix;

    // Extract derivative matrix of current coefficient matrix
    Eigen::MatrixXd currentCoefficientMatrix = spacecraftPointer_->getCurrentCoefficientMatrix();
    double v_r = velocity_e1;
    double dv_r = currentCoefficientMatrix( 1, 8 ); // 0th row is current state, 1st row is first derivative of state
    double ddv_r = currentCoefficientMatrix( 2, 8 );
    double dddv_r = currentCoefficientMatrix( 3, 8 );
    double v_tr = velocity_e2;
    double dv_tr = currentCoefficientMatrix( 1, 9 ); // 0th row is current state, 1st row is first derivative of state
    double ddv_tr = currentCoefficientMatrix( 2, 9 );
    double dddv_tr = currentCoefficientMatrix( 3, 9 );

    // Analytic calculation of the first three derivatives of the flight path angle gamma
    double gammaDot = ( dv_r /v_tr - (v_r*dv_tr)/(v_tr * v_tr))/((v_r * v_r)/(v_tr * v_tr) + 1.0 );
    double gammaDDot = (ddv_r/v_tr - (2.0*dv_r*dv_tr)/(v_tr * v_tr)
                        + (2.0*v_r*dv_tr * dv_tr)/(v_tr * v_tr * v_tr) - (v_r*ddv_tr)/(v_tr * v_tr))/((v_r * v_r)/(v_tr * v_tr) + 1.0)
            - ((dv_r/v_tr - (v_r*dv_tr)/(v_tr * v_tr))
               *((2.0*v_r*dv_r)/(v_tr * v_tr) - (2.0*(v_r * v_r)*dv_tr)/(v_tr * v_tr * v_tr)))/
            std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 2.0 );
    double gammaDDDot = (2.0*(dv_r/v_tr
                            - (v_r*dv_tr)/(v_tr * v_tr))*std::pow( ((2.0*v_r*dv_r)/(v_tr * v_tr)
                                                                    - (2.0*(v_r * v_r)*dv_tr)/(v_tr * v_tr * v_tr)), 2.0) )/
            std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 3.0)
            - (2.0*((2.0*v_r*dv_r)/(v_tr * v_tr) - (2.0*(v_r * v_r)*dv_tr)/(v_tr * v_tr * v_tr))*(ddv_r/v_tr
                                                                               - (2.0*dv_r*dv_tr)/(v_tr * v_tr)
                                                                               + (2.0*v_r*dv_tr * dv_tr)/(v_tr * v_tr * v_tr)
                                                                               - (v_r*ddv_tr)/(v_tr * v_tr)))/
            std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 2.0 )
            - ((dv_r/v_tr - (v_r*dv_tr)/(v_tr * v_tr))*((2.0*dv_r * dv_r)/(v_tr * v_tr)
                                                    + (6.0*(v_r * v_r)*dv_tr * dv_tr)/(v_tr * v_tr * v_tr * v_tr)
                                                    - (2.0*(v_r * v_r)*ddv_tr)/(v_tr * v_tr * v_tr)
                                                    + (2.0*v_r*ddv_r)/(v_tr * v_tr)
                                                    - (8.0*v_r*dv_r*dv_tr)/(v_tr * v_tr * v_tr)))
            /std::pow( ((v_r * v_r)/(v_tr * v_tr) + 1.0), 2.0 ) - ((v_r*dddv_tr)/(v_tr * v_tr) - dddv_r/v_tr
                                           - (6.0*dv_r*dv_tr * dv_tr)/(v_tr * v_tr * v_tr) + (3.0*dv_r*ddv_tr)/(v_tr * v_tr)
                                           + (3.0*dv_tr*ddv_r)/(v_tr * v_tr) + (6.0*v_r*dv_tr * dv_tr * dv_tr)/(v_tr * v_tr * v_tr * v_tr)
                                           - (6.0*v_r*dv_tr*ddv_tr)/(v_tr * v_tr * v_tr))/((v_r * v_r)/(v_tr * v_tr) + 1.0);

    Eigen::Matrix< double, 3, 3> transformationMatrixVftoUSMframe;
    transformationMatrixVftoUSMframe <<
            std::sin( gamma ), -std::cos( gamma ), 0.0,
            std::cos( gamma ), std::sin( gamma ), 0.0,
            0.0, 0.0, 1.0;
    Eigen::Matrix< double, 3, 3> firstDerivativeOfTransformationMatrixVfToUSMframe;
    firstDerivativeOfTransformationMatrixVfToUSMframe <<
            gammaDot * std::cos( gamma ), gammaDot * std::sin( gamma ), 0.0,
            -gammaDot * sin( gamma ), gammaDot * std::cos( gamma ), 0.0,
            0.0, 0.0, 0.0;
    Eigen::Matrix< double, 3, 3> secondDerivativeOfTransformationMatrixVfToUSMframe;
    secondDerivativeOfTransformationMatrixVfToUSMframe <<
            gammaDDot * std::cos( gamma ) - gammaDot * gammaDot * std::sin( gamma ),
            gammaDDot * std::sin( gamma ) + gammaDot * gammaDot * std::cos( gamma ), 0.0,
            -gammaDDot * std::sin( gamma ) - gammaDot * gammaDot * std::cos( gamma ),
            gammaDDot * std::cos( gamma ) - gammaDot * gammaDot * std::sin( gamma ), 0.0,
            0.0, 0.0, 0.0;
    Eigen::Matrix< double, 3, 3> thirdDerivativeOfTransformationMatrixVfToUSMframe;
    thirdDerivativeOfTransformationMatrixVfToUSMframe <<
            gammaDDDot * std::cos( gamma )
            - 3.0 * gammaDDot * gammaDot * std::sin( gamma )
            - gammaDot * gammaDot * gammaDot * std::cos( gamma ),
            gammaDDDot * std::sin( gamma )
            + 3.0 * gammaDDot * gammaDot * std::cos( gamma )
            - gammaDot * gammaDot * gammaDot * std::sin( gamma ), 0.0,
            - gammaDDDot * std::sin( gamma )
            - 3.0 * gammaDDot * gammaDot * std::cos( gamma )
            + gammaDot * gammaDot * gammaDot * std::sin( gamma ),
            gammaDDDot * std::cos( gamma )
            - 3.0 * gammaDDot * gammaDot * std::sin( gamma )
            - gammaDot * gammaDot * gammaDot * std::cos( gamma ), 0.0,
            0.0, 0.0, 0.0;


    // First derivative
    currentDerivativesOfThrustInUSMframe.row( 0 ) =
            ( firstDerivativeOfTransformationMatrixVfToUSMframe * currentThrustForceInVelocityFrame.adjoint( )
              + transformationMatrixVftoUSMframe * currentDerivativesOfThrustInVelocityFrame.row( 0 ).adjoint( ) ).adjoint( );
    currentDerivativesOfThrustInUSMframe.row( 1 ) =
            ( secondDerivativeOfTransformationMatrixVfToUSMframe * currentThrustForceInVelocityFrame.adjoint( )
              + 2.0 * firstDerivativeOfTransformationMatrixVfToUSMframe *
                      currentDerivativesOfThrustInVelocityFrame.row( 0 ).adjoint( )
              + transformationMatrixVftoUSMframe * currentDerivativesOfThrustInVelocityFrame.row( 1 ).adjoint( ) ).adjoint( );
    currentDerivativesOfThrustInUSMframe.row( 2 ) =
            ( thirdDerivativeOfTransformationMatrixVfToUSMframe * currentThrustForceInVelocityFrame.adjoint( )
              + 3.0 * secondDerivativeOfTransformationMatrixVfToUSMframe *
                      currentDerivativesOfThrustInVelocityFrame.row( 0 ).adjoint( )
              + 3.0 * firstDerivativeOfTransformationMatrixVfToUSMframe *
                      currentDerivativesOfThrustInVelocityFrame.row( 1 ).adjoint( )
              + transformationMatrixVftoUSMframe * currentDerivativesOfThrustInVelocityFrame.row( 2 ).adjoint( ) ).adjoint( );

    //??
    // std::cout << "label 6" << std::endl;

    // Calculate derivatives of current thrust acceleration in USM frame components
    Eigen::VectorXd currentMassDerivatives( 3 ); // first three derivatives of mass
    currentMassDerivatives( 0 ) = - currentResultantThrustForce
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );
    currentMassDerivatives( 1 ) = - currentDerivativesOfResultantThrustForce( 0, 0 )
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );
    currentMassDerivatives( 2 ) = - currentDerivativesOfResultantThrustForce( 1, 0 )
            / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse() );

    //??
    // // std::cout << "label 7" << std::endl;

    Eigen::MatrixXd currentDerivativesOfThrustAccelerationInUSMframe( std::max( order, 3 ), 3 ); // At least 3 derivatives are calculated
    currentDerivativesOfThrustAccelerationInUSMframe.fill( 0.0 );

    //??
    // // std::cout << "label 8" << std::endl;

    for ( int j = 0; j < 3; j++ )
    {
        // First derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 0, j ) =
                ( currentDerivativesOfThrustInUSMframe( 0, j ) * currentMass
                  - currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass );
        // Second derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 1, j ) =
                ( currentDerivativesOfThrustInUSMframe( 1, j ) * currentMass * currentMass
                  - currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMass
                  - 2.0 * currentMass * currentMassDerivatives( 0 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                  + 2.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass * currentMass );
        // Third derivative
        currentDerivativesOfThrustAccelerationInUSMframe( 2, j ) =
                ( ( currentDerivativesOfThrustInUSMframe( 2, j ) * currentMass * currentMass
                    - currentMassDerivatives( 2 ) * currentThrustForceInUSMframe( j ) * currentMass
                    - currentMassDerivatives( 1 ) * currentDerivativesOfThrustInUSMframe( 0, j ) * currentMass
                    - currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMassDerivatives( 0 )
                    - 2.0 * currentMass * currentMassDerivatives( 1 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                    + 4.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) )
                  * currentMass
                  - 3.0 * currentDerivativesOfThrustInUSMframe( 1, j ) * currentMass * currentMass
                  + 3.0 * currentMassDerivatives( 1 ) * currentThrustForceInUSMframe( j ) * currentMass
                  + 6.0 * currentMass * currentMassDerivatives( 0 ) * currentDerivativesOfThrustInUSMframe( 0, j )
                  - 6.0 * currentMassDerivatives( 0 ) * currentMassDerivatives( 0 ) * currentThrustForceInUSMframe( j ) )
                / ( currentMass * currentMass * currentMass * currentMass );
    }


    // Save current derivatives of thrust acceleration in USM frame components in spacecraft pointer
    spacecraftPointer_->setCurrentAccelerationDerivative( currentDerivativesOfThrustAccelerationInUSMframe );

}
