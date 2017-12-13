#include "constants.h"
#include "spacecraft.h"
#include "problemRecurrenceRelations.h"
#include "recurrenceRelations.h"
#include "spacecraft.h"

#include "recurrenceRelations.h"

#include "Eigen/Core"

#include <boost/math/special_functions/factorials.hpp>
#include <boost/make_shared.hpp>


ProblemRecurrenceRelations::ProblemRecurrenceRelations(
        int order,
        SpacecraftPointer spacecraftPointer,
        ConstantsPointer constantsPointer ):
        order_( order ),
        spacecraftPointer_( spacecraftPointer ),
        constantsPointer_( constantsPointer ) {}

void ProblemRecurrenceRelations::computeNextState(
        double currentStepSize )
{
    // This function carries out the state advancement in a central
    // body gravitational field with arbitrary thrust input

    Vector11d currentState = spacecraftPointer_->getCurrentState();
    // Define number of state variables
    int numberOfStateVariables = currentState.rows(); // Assuming that currentState is a columnvector

    // // Initialize the state
    double x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x7_0, x8_0, x9_0,  x10_0, x11_0;

    x1_0 = currentState( 0 );
    x2_0 = currentState( 1 );
    x3_0 = currentState( 2 );
    x4_0 = currentState( 3 );
    x5_0 = currentState( 4 );
    x6_0 = currentState( 5 );
    x7_0 = currentState( 6 );
    x8_0 = currentState( 7 );
    x9_0 = currentState( 8 );
    x10_0 = currentState( 9 );
    x11_0 = currentState( 10 );

    // // Declare some vectors
    // State
    Eigen::VectorXd X1( order_ + 1 ), // Use dynamic matrix instead of fixed Eigen::Matrix < double, order_ + 1, 1 >
                    X2( order_ + 1 ),
                    X3( order_ + 1 ),
                    X4( order_ + 1 ),
                    X5( order_ + 1 ),
                    X6( order_ + 1 ),
                    X7( order_ + 1 ),
                    X8( order_ + 1 ),
                    X9( order_ + 1 ),
                    X10( order_ + 1 ),
                    X11( order_ + 1 );
    // Time
    Eigen::VectorXd T_vector( order_ + 1 );
    // Auxiliary
    Eigen::VectorXd W1( order_ + 1 ),
                    W2_1( order_ + 1 ),
                    W2_2( order_ + 1 ),
                    W3_1( order_ + 1 ),
                    W3_2( order_ + 1 ),
                    W4_1( order_ + 1 ),
                    W4_2( order_ + 1 ),
                    W5_1( order_ + 1 ),
                    W5_2( order_ + 1 ),
                    W6_1( order_ + 1 ),
                    W6_2( order_ + 1 ),
                    W7_1( order_ + 1 ),
                    W7_2( order_ + 1 ),
                    W9_1( order_ + 1 ),
                    W9_2( order_ + 1 ),
                    W10( order_ + 1 ),
                    W11_1( order_ + 1 ),
                    W11_2( order_ + 1 );

    // // Time derivatives
    Eigen::VectorXd U1( order_ + 1 ),
                    U2( order_ + 1 ),
                    U3( order_ + 1 ),
                    U4( order_ + 1 ),
                    U5( order_ + 1 ),
                    U6( order_ + 1 ),
                    U7( order_ + 1 ),
                    U8( order_ + 1 ),
                    U9( order_ + 1 ),
                    U10( order_ + 1 ),
                    U11( order_ + 1 );

    // // Auxiliary V variables
    Eigen::VectorXd V1_1( order_ + 1 ), // q3^2
                    V1_2( order_ + 1 ), // q4^2
                    V1( order_ + 1 ), // q3^2 + q4^2
                    V2( order_ + 1 ), // q3 * q4
                    V3( order_ + 1 ), // q4^2 - q3^2
                    V4_1( order_ + 1 ), // q1 * q3
                    V4_2( order_ + 1 ), // q2 * q4
                    V4( order_ + 1 ),
                    V5_1( order_ + 1 ), // V2/V1
                    V5( order_ + 1 ), // sin(lambda)
                    V6( order_ + 1 ), // V3/V1 = cos(lambda)
                    V7( order_ + 1 ), // gamma
                    V8( order_ + 1 ), // C/V_e2 = p
                    V9( order_ + 1 ), // gamma/Ve2
                    V10( order_ + 1 ), // Rf1 * gamma/ve2
                    V11( order_ + 1 ), // Rf2 * gamma/ve2
                    V12( order_ + 1 ), // ae1
                    V13( order_ + 1 ), // ae2
                    V14( order_ + 1 ), // ae3
                    V15( order_ + 1 ), // w1
                    V16( order_ + 1 ), // ae2(1+p)
                    V17( order_ + 1 ), // ae3 * gamma * Rf1/ve2
                    V18( order_ + 1 ); // ae3 * gamma * Rf2/ve2

    // // Giving the starting conditions
    X1( 0 ) = x1_0;
    X2( 0 ) = x2_0;
    X3( 0 ) = x3_0;
    X4( 0 ) = x4_0;
    X5( 0 ) = x5_0;
    X6( 0 ) = x6_0;
    X7( 0 ) = x7_0;
    X8( 0 ) = x8_0;
    X9( 0 ) = x9_0;
    X10( 0 ) = x10_0;
    X11( 0 ) = x11_0;

    // Specify current mass
    double m = x8_0;

    // h^k = 1 when k = 0
    T_vector( 0 ) = 1.0;


    // // Auxiliary V variables
    V1_1( 0 ) = X6( 0 ) * X6( 0 ); // q3^2
    V1_2( 0 ) = X7( 0 ) * X7( 0 ); // q4^2
    V1( 0 ) = V1_1( 0 ) + V1_2( 0 );  // q3^2 + q4^2

    V2( 0 ) = X6( 0 ) * X7( 0 ); // q3 * q4

    V3( 0 ) = V1_2( 0 ) - V1_1( 0 ); // q4^2 - q3^2

    V4_1( 0 ) = X4( 0 ) * X6( 0 ); // q1 * q3
    V4_2( 0 ) = X5( 0 ) * X7( 0 );// q2 * q4;
    V4( 0 ) = V4_1( 0 ) - V4_2( 0 );


    V5_1( 0 ) = V2( 0 ) / V1( 0 ); // V2/V1
    V5( 0 ) = 2 * V5_1( 0 );  // sin(lambda)

    V6( 0 ) = V3( 0 ) / V1( 0 );  // V3/V1 = cos(lambda)

    V7( 0 ) = V4( 0 ) / V1( 0 ); // gamma

    V8( 0 ) = X1( 0 ) / X10( 0 ); // C/V_e2 = p

    V9( 0 ) = V7( 0 ) / X10( 0 ); // gamma/Ve2

    V10( 0 ) = X2( 0 ) * V9( 0 ); // Rf1 * gamma/ve2

    V11( 0 ) = X3( 0 ) * V9( 0 ); // Rf2 * gamma/ve2

    // // New ones
    // get current acceleration from spacecraft class
    Eigen::Matrix< double, 1, 3 > acc_curr;
    acc_curr = spacecraftPointer_->getCurrentAcceleration();

    //??
//    std::cout << "acc_curr = \n" << acc_curr << std::endl;

    V12( 0 ) = acc_curr( 0 ); // ae1

    V13( 0 ) = acc_curr( 1 ); // ae2

    V14( 0 ) = acc_curr( 2 ); // ae3


    V15( 0 ) = V14( 0 ) / X10( 0 ); // omega 1

    V16( 0 ) = V13( 0 ) + V13( 0 ) * V8( 0 ); // ae2(1+p)

    V17( 0 ) = V14( 0 ) * V10( 0 ); // ae3*gamma*Rf1/ve2

    V18( 0 ) = V14( 0 ) * V11( 0 ); // ae3*gamma*Rf2/ve2

    // // Time Derivatives and Auxilary Functions

    W1( 0 ) = V8( 0 ) * V13( 0 );
    U1( 0 ) = -W1( 0 );


    W2_1( 0 ) = V12( 0 ) * V6( 0 );
    W2_2( 0 ) = V16( 0 ) * V5( 0 );
    U2( 0 ) = W2_1( 0 ) - W2_2( 0 ) - V18( 0 );


    W3_1( 0 ) = V12( 0 ) * V5( 0 );
    W3_2( 0 ) = V16( 0 ) * V6( 0 );
    U3( 0 ) = W3_1( 0 ) + W3_2( 0 ) + V17( 0 );

    W4_1( 0 ) = X11( 0 ) * X5( 0 );
    W4_2( 0 ) = V15( 0 ) * X7( 0 );
    U4( 0 ) = 0.5 * W4_1( 0 ) + 0.5 * W4_2( 0 );

    W5_1( 0 ) = X11( 0 ) * X4( 0 );
    W5_2( 0 ) = V15( 0 ) * X6( 0 );
    U5( 0 ) = -0.5 * W5_1( 0 ) + 0.5 * W5_2( 0 );

    W6_1( 0 ) = X11( 0 ) * X7( 0 );
    W6_2( 0 ) = V15( 0 ) * X5( 0 );
    U6( 0 ) = 0.5 * W6_1( 0 ) - 0.5 * W6_2( 0 );

    W7_1( 0 ) = X11( 0 ) * X6( 0 );
    W7_2( 0 ) = V15( 0 ) * X4( 0 );
    U7( 0 ) = -0.5 * W7_1( 0 ) - 0.5 * W7_2( 0 );

    Eigen::Matrix< double, 1, 3 > currentThrust = acc_curr * m;
    U8( 0 ) = -currentThrust.norm() /
            ( constantsPointer_->standardGravity_
              * spacecraftPointer_->getSpecificImpulse() ); // Current thrust force = constant acc * current mass

    //??
//    std::cout << "currentThrust = " << currentThrust << std::endl;

    W9_1( 0 ) = X1( 0 ) * X11( 0 );
    W9_2( 0 ) = X10( 0 ) * X11( 0 );
    U9( 0 ) = V12( 0 ) - W9_1( 0 ) + W9_2( 0 );

    W10( 0 ) = X9( 0 ) * X11( 0 );
    U10( 0 ) = V13( 0 ) - W10( 0 );

    W11_1( 0 ) = X11( 0 ) * U1( 0 ) / X1( 0 );
    W11_2( 0 ) = X11( 0 ) * U10( 0 ) / X10( 0 );
    U11( 0 ) = W11_1( 0 ) + 2 * W11_2( 0 );

    Eigen::MatrixXd dacc_curr( std::max( order_, 3 ), 3 );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > dT_res;

    // // This is getting the Taylor series
    for ( int k = 1; k <= order_; k++ )
    {
        // Make object of class ProblemRecurrenceRelations and set order
        RecurrenceRelationsPointer recurrenceRelationsPointer =
                boost::make_shared< RecurrenceRelations >( k );

        V1_1( k ) = recurrenceRelationsPointer->multiplyRecursive( X6, U6, X6, U6 ); //X6(1)*X6(1);
        V1_2( k ) = recurrenceRelationsPointer->multiplyRecursive( X7, U7, X7, U7 ); //X7(1)*X7(1);
        V1( k ) = V1_1( k ) + V1_2( k );

        V2( k ) = recurrenceRelationsPointer->multiplyRecursive( X6, U6, X7, U7 ); //X6(1)*X7(1);

        V3( k ) = V1_2( k ) - V1_1( k ); //q4^2 - q3^2

        V4_1( k ) = recurrenceRelationsPointer->multiplyRecursive( X4, U4, X6, U6 ); //X4(1)*X6(1); //q1*q3
        V4_2( k ) = recurrenceRelationsPointer->multiplyRecursive( X5, U5, X7, U7 ); //X5(1)*X7(1);//q2*q4;
        V4( k ) = V4_1( k ) - V4_2( k );

        V5_1( k ) = recurrenceRelationsPointer->V1overV2( V5_1, V2, V1 ); // V2(1)/V1(1); //V2/V1
        V5( k ) = 2 * V5_1( k ); //sin(lambda)

        V6( k ) = recurrenceRelationsPointer->V1overV2( V6, V3, V1 ); // V3(1)/V1(1);  //V3/V1 = cos(lambda)

        V7( k ) = recurrenceRelationsPointer->V1overV2( V7, V4, V1 ); // V4(1)/V1(1); //gamma

        V8( k ) = recurrenceRelationsPointer->divideRecursive( V8, U1, U10, X10 ); //X1(1)/X10(1); //C/V_e2 = p

        V9( k ) = recurrenceRelationsPointer->Um_over_Xn( V9, V7, U10, X10 ); //V7(1)/X10(1); //gamma/Ve2

        V10( k ) = recurrenceRelationsPointer->X1xV2_Recursive( X2, U2, V9 ); //X2(1)*V9(1); //Rf1*gamma/ve2

        V11( k ) = recurrenceRelationsPointer->X1xV2_Recursive( X3, U3, V9 ); //X3(1)*V9(1); //Rf2*gamma/ve2


        dacc_curr = spacecraftPointer_->getCurrentAccelerationDerivative();


        V12( k ) = dacc_curr( k - 1, 0 ); //ae1 // no use of k because k=1 corresponds to second derivative already

        V13( k ) = dacc_curr( k - 1, 1 ); //ae2

        V14( k ) = dacc_curr( k - 1, 2 ); //ae3


        V15( k ) = recurrenceRelationsPointer->Um_over_Xn( V15, V14, U10, X10 ); //V14(1)/X10(1); //omega 1

        V16( k ) = V13( k ) + recurrenceRelationsPointer->V1xV2( V13, V8 ); //V13(1)*V8(1);//ae2(1+p)

        V17( k ) = recurrenceRelationsPointer->V1xV2( V14, V10 ); //V14(1)*V10(1); //ae3*gamma*Rf1/ve2

        V18( k ) = recurrenceRelationsPointer->V1xV2( V14, V11 ); //V14(1)*V11(1); //ae3*gamma*Rf2/ve2


        // Compute the auxiliary functions first
        W1( k ) = recurrenceRelationsPointer->V1xV2( V8, V13 ); //V14(1)*V11(1); //ae3*gamma*Rf2/ve2
        U1( k ) = -W1( k );

        W2_1( k ) = recurrenceRelationsPointer->V1xV2( V12, V6 ); //V12(1)*V6(1);
        W2_2( k ) = recurrenceRelationsPointer->V1xV2( V16, V5 ); //V16(1)*V5(1);
        U2( k ) = W2_1( k ) - W2_2( k ) - V18( k );

        W3_1( k ) = recurrenceRelationsPointer->V1xV2( V12, V5 ); //V12(1)*V5(1);
        W3_2( k ) = recurrenceRelationsPointer->V1xV2( V16, V6 );//V16(1)*V6(1); //V16(1)*V6(1);
        U3( k ) = W3_1( k ) + W3_2( k ) + V17( k );

        W4_1( k ) = recurrenceRelationsPointer->multiplyRecursive( X11, U11, X5, U5 );//X11(1)*X5(1);
        W4_2( k ) = recurrenceRelationsPointer->X1xV2_Recursive( X7, U7, V15 );//V15(1)*X7(1);
        U4( k ) = 0.5 * W4_1( k ) + 0.5 * W4_2( k );

        W5_1( k ) = recurrenceRelationsPointer->multiplyRecursive( X11, U11, X4, U4 );//X11(1)*X4(1);
        W5_2( k ) = recurrenceRelationsPointer->X1xV2_Recursive( X6, U6, V15 );//V19(1)*X6(1);
        U5( k ) = -0.5 * W5_1( k ) + 0.5 * W5_2( k );

        W6_1( k ) = recurrenceRelationsPointer->multiplyRecursive( X11, U11, X7, U7 );//X11(1)*X7(1);
        W6_2( k ) = recurrenceRelationsPointer->X1xV2_Recursive( X5, U5, V15 );//V15(1)*X5(1);
        U6( k ) = 0.5 * W6_1( k ) - 0.5 * W6_2( k );

        W7_1( k ) = recurrenceRelationsPointer->multiplyRecursive( X11, U11, X6, U6 );//X11(1)*X6(1);
        W7_2( k ) = recurrenceRelationsPointer->X1xV2_Recursive( X4, U4, V15 );//V15(1)*X4(1);
        U7( k ) = -0.5 * W7_1( k ) - 0.5 * W7_2( k );

        /// Equation added by Michael Van den Broeck
        // Get current resultant low-thrust force derivatives
        dT_res = spacecraftPointer_->getCurrentResultantThrustForceDerivatives();

        U8( k ) = - dT_res( k - 1 ) / ( constantsPointer_->standardGravity_ * spacecraftPointer_->getSpecificImpulse()
                                    * boost::math::factorial<double>( k ) ); // no use of  k because k=1 corresponds to second derivative already
        ///

        W9_1( k ) = recurrenceRelationsPointer->multiplyRecursive( X1, U1, X11, U11 );//X1(1)*X11(1);
        W9_2( k ) = recurrenceRelationsPointer->multiplyRecursive( X10, U10, X11, U11 );//X10(1)*X11(1);
        U9( k ) = V12( k ) - W9_1( k ) + W9_2( k );

        W10( k ) = recurrenceRelationsPointer->multiplyRecursive( X9, U9, X11, U11 );//X9(1)*X11(1);
        U10( k ) = V13( k ) - W10( k );

        W11_1( k ) = recurrenceRelationsPointer->X1xU2overX2_Recursive( W11_1, X11, U11, X1, U1 );//X11(1)*U1(1)/X1(1);
        W11_2( k ) = recurrenceRelationsPointer->X1xU2overX2_Recursive( W11_2, X11, U11, X10, U10 );//X11(1)*U10(1)/X10(1);
        U11( k ) = W11_1( k ) + 2 * W11_2( k );






        //// ALL THE X ////
        X1( k ) = U1( k - 1) / k;  // UpdateX_from_U(U1,k); // dim( X1 ) = (order+1) x 1
        X2( k ) = U2( k - 1) / k;  // UpdateX_from_U(U2,k);
        X3( k ) = U3( k - 1) / k;  // UpdateX_from_U(U3,k);
        X4( k ) = U4( k - 1) / k;  // UpdateX_from_U(U4,k);
        X5( k ) = U5( k - 1) / k;  // UpdateX_from_U(U5,k);
        X6( k ) = U6( k - 1) / k;  // UpdateX_from_U(U6,k);
        X7( k ) = U7( k - 1) / k;  // UpdateX_from_U(U7,k);
        X8( k ) = U8( k - 1) / k;  // UpdateX_from_U(U8,k);
        X9( k ) = U9( k - 1) / k;  // UpdateX_from_U(U9,k);
        X10( k ) = U10( k - 1) / k;  // UpdateX_from_U(U10,k);
        X11( k ) = U11( k - 1) / k;  // UpdateX_from_U(U11,k);


        T_vector( k ) = std::pow( currentStepSize, k ); // dim(T_vector) = (order+1) x 1

    }

    // Declare variables
    Eigen::MatrixXd x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9, x_10, x_11;

    x_1 = X1.adjoint() * T_vector; // scalar
    x_2 = X2.adjoint() * T_vector; // scalar
    x_3 = X3.adjoint() * T_vector; // scalar
    x_4 = X4.adjoint() * T_vector; // scalar
    x_5 = X5.adjoint() * T_vector; // scalar
    x_6 = X6.adjoint() * T_vector; // scalar
    x_7 = X7.adjoint() * T_vector; // scalar
    x_8 = X8.adjoint() * T_vector; // scalar
    x_9 = X9.adjoint() * T_vector; // scalar
    x_10 = X10.adjoint() * T_vector; // scalar
    x_11 = X11.adjoint() * T_vector; // scalar

    // Return the next state
    Eigen::MatrixXd State_out( numberOfStateVariables, 1 );
    State_out << x_1( 0 ), x_2( 0 ), x_3( 0 ), x_4( 0 ), x_5( 0 ),
            x_6( 0 ), x_7( 0 ), x_8( 0 ), x_9( 0 ), x_10( 0 ), x_11( 0 );

    // Store current state vector in spacecraft pointer
    spacecraftPointer_->setCurrentState( State_out );

    Eigen::MatrixXd X_matrix( order_ + 1, numberOfStateVariables);
    X_matrix << X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11;

    // Store current coefficient matrix of derivatives in spacecraft pointer
    spacecraftPointer_->setCurrentCoefficientMatrix( X_matrix );

}

ProblemRecurrenceRelations::~ProblemRecurrenceRelations()
{

}

