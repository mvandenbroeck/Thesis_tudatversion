#include "recurrenceRelations.h"
#include "Eigen/Core"

double RecurrenceRelations::multiplyRecursive(
        const Eigen::VectorXd &X_left,                  // dim: (order+1) x 1
        const Eigen::VectorXd &U_left,
        const Eigen::VectorXd &X_right,
        const Eigen::VectorXd &U_right )
{
    // Find the product
    // W_P = X_left * X_right.

    int k = orderInLoopOfTaylorSeries_;

    double W_P_k = X_left( 0 ) * U_right( k - 1 ) / k +
            X_right( 0 ) * U_left( k - 1 ) / k;

    double Sum_positive, Sum_negative;
    Sum_positive = 0.0;
    Sum_negative = 0.0;

    int k_minus_j;
    double Value;

    if ( k > 1 )
    {
        for ( int j = 1; j < k; j++ )
        {
            k_minus_j = k - j;
            Value = ( U_left( j - 1 ) / j ) * ( U_right( k_minus_j - 1 ) / ( k - j ) );

            if ( Value > 0)
            {
                Sum_negative = Sum_negative + Value;
            }
            else
            {
                Sum_positive = Sum_positive + Value;
            }
        }
    }

    double Sum;
    Sum = Sum_positive + Sum_negative; // To avoid loss of significant digits

    W_P_k = W_P_k + Sum;
    return W_P_k;
}


double RecurrenceRelations::V1overV2(
        const Eigen::VectorXd &W,               // dim: (order+1) x 1
        const Eigen::VectorXd &V1,
        const Eigen::VectorXd &V2 )
{
    // W = V1/V2

    int k = orderInLoopOfTaylorSeries_;
    double W_k = V1( k );

    double Sum_positive, Sum_negative;
    Sum_positive = 0.0;
    Sum_negative = 0.0;

    double Value;

    if ( k > 0 )
    {
        for ( int j = 1; j <= k; j++ )
        {
            Value = V2( j ) * W( k - j );

            if ( Value < 0 )
            {
                Sum_negative = Sum_negative + Value;
            }
            else
            {
                Sum_positive = Sum_positive + Value;
            }
        }
    }

    double Sum = Sum_positive + Sum_negative; // To avoid loss of significant digits

    W_k = ( W_k - Sum ) / V2( 0 );

    return W_k;

}


double RecurrenceRelations::divideRecursive(
        const Eigen::VectorXd &W,               // dim: (order+1) x 1
        const Eigen::VectorXd &U_num,
        const Eigen::VectorXd &U_denom,
        const Eigen::VectorXd &X_denom )
{
    // w = x_num / x_denom

    int k = orderInLoopOfTaylorSeries_;

    double W_k = ( U_num( k - 1 ) / k );

    double Sum_positive, Sum_negative;
    Sum_positive = 0.0;
    Sum_negative = 0.0;

    double Value;

    if ( k > 0 )
    {
        for ( int j = 1; j <= k; j++ )
        {
            Value = ( U_denom( j - 1 ) / j ) * W( k - j );

            if ( Value < 0 )
            {
                Sum_negative = Sum_negative + Value;
            }
            else
            {
                Sum_positive = Sum_positive + Value;
            }
        }
    }

    double Sum = Sum_positive + Sum_negative; // To avoid loss of significant digits

    W_k = ( W_k - Sum ) / X_denom( 0 );

    return W_k;

}

double RecurrenceRelations::Um_over_Xn(
        const Eigen::VectorXd &W,          // dim: (order+1) x 1
        const Eigen::VectorXd &U_num,
        const Eigen::VectorXd &U_denom,
        const Eigen::VectorXd &X_denom )
{
    // w = u_num / x_denom

    int k = orderInLoopOfTaylorSeries_;

    double W_k = U_num( k );

    double Sum_positive, Sum_negative;
    Sum_positive = 0.0;
    Sum_negative = 0.0;

    double Value;

    if ( k > 0 )
    {
        for ( int j = 1; j <= k; j++ )
        {
            Value = ( U_denom( j - 1 ) / j ) * W( k - j );

            if ( Value < 0 )
            {
                Sum_negative = Sum_negative + Value;
            }
            else
            {
                Sum_positive = Sum_positive + Value;
            }
        }
    }

    double Sum = Sum_positive + Sum_negative; // To avoid loss of significant digits

    W_k = ( W_k - Sum ) / X_denom( 0 );

    return W_k;
}


double RecurrenceRelations::X1xV2_Recursive(
        const Eigen::VectorXd &X_L,          // dim: (order+1) x 1
        const Eigen::VectorXd &U_L,
        const Eigen::VectorXd &V_R )
{
    // w = x_L * u_R

    int k = orderInLoopOfTaylorSeries_;

    double Part1 = X_L( 0 ) * V_R( k );

    double Part2_positive = 0.0;
    double Part2_negative = 0.0;

    int k_minus_j;
    double Value;

    for ( int j = 1; j <= k; j++ )
    {
        k_minus_j = k - j;
        Value = ( U_L( j - 1 )/j ) * V_R( k_minus_j );

        if ( Value < 0 )
        {
            Part2_negative = Part2_negative + Value;
        }
        else
        {
            Part2_positive = Part2_positive + Value;
        }
    }

    double Part2 = Part2_positive + Part2_negative;

    double W_k = Part1 + Part2;

    return W_k;

}


double RecurrenceRelations::V1xV2(
        const Eigen::VectorXd &V1,          // dim: (order+1) x 1
        const Eigen::VectorXd &V2 )
{
    int k = orderInLoopOfTaylorSeries_;

    double Sum_positive, Sum_negative;
    Sum_positive = 0.0;
    Sum_negative = 0.0;

    int k_minus_j;
    double Value;

    for ( int j = 0; j <= k; j++ )
    {
        k_minus_j = k - j;

        Value = V1( j ) * V2( k_minus_j );

        if ( Value < 0 )
        {
            Sum_negative = Sum_negative + Value;
        }
        else
        {
            Sum_positive = Sum_positive + Value;
        }
    }

    double Sum = Sum_positive + Sum_negative; // To avoid loss of significant digits

    double W_P_k = Sum;

    return W_P_k;
}


double RecurrenceRelations::X1xU2overX2_Recursive(const Eigen::VectorXd &W,          // dim: (order+1) x 1
        const Eigen::VectorXd &X_L,
        const Eigen::VectorXd &U_L,
        const Eigen::VectorXd &X_R,
        const Eigen::VectorXd &U_R)
{
    // w = x_L * u_R / x_R

    int k = orderInLoopOfTaylorSeries_;

    double Part1 = X_L( 0 ) * U_R( k );

    double Part2_positive = 0.0;
    double Part2_negative = 0.0;

    double Part3_positive = 0.0;
    double Part3_negative = 0.0;

    int k_minus_j;
    double Value;

    for ( int j = 1; j <= k; j++ )
    {
        k_minus_j = k - j;
        Value = ( U_L( j - 1 ) /j ) * U_R( k_minus_j );

        if ( Value < 0 )
        {
            Part2_negative = Part2_negative + Value;
        }
        else
        {
            Part2_positive = Part2_positive + Value;
        }

        Value = ( U_R( j - 1 ) / j) * W( k_minus_j );

        if ( Value < 0 )
        {
            Part3_negative = Part3_negative + Value;
        }
        else
        {
            Part3_positive = Part3_positive + Value;
        }
    }

    double Part2 = Part2_positive + Part2_negative;

    double Part3 = Part3_positive + Part3_negative;

    double W_k = ( 1 / X_R( 0 ) ) * ( Part1 + Part2 - Part3 );

    return W_k;

}
