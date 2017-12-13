#include "discreteForceModel.h"

#include <Eigen/Core>

#include <boost/make_shared.hpp>

#include "spacecraft.h"
#include "integrationSettings.h"
#include "constants.h"
#include "lagrangeInterpolator.h"
#include "frameTransformation.h"
#include "tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"

DiscreteForceModel::DiscreteForceModel(
        SpacecraftPointer spacecraftPointer,
        IntegrationSettingsPointer integrationSettingsPointer,
        ConstantsPointer constantsPointer ):
        spacecraftPointer_( spacecraftPointer ),
        integrationSettingsPointer_( integrationSettingsPointer ),
        constantsPointer_( constantsPointer ) {}

void DiscreteForceModel::updateCurrentForcesAndAccelerationsForTSI(
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

DiscreteForceModel::~DiscreteForceModel()
{

}

