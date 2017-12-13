#include <Eigen/Core>
#include "lagrangeInterpolator.h"




LagrangeInterpolator::LagrangeInterpolator(
            Eigen::Matrix< double, Eigen::Dynamic, 2 > pointsToBeInterpolated ):
            pointsToBeInterpolated_( pointsToBeInterpolated ) {}



Eigen::Matrix< double, Eigen::Dynamic, 1> LagrangeInterpolator::computeCoefficientsOfInterpolatingPolynomial( )
{
/*  This function creates a Lagrange polynomial that interpolates the vector
    of input data. It builds a polynomial of order n, in which n is the
    number of data points in the input vector. The output of the function is
    the vector containing the coefficients of the Lagrange polynomial in
    decreasing order of x.

    Input:
      PointsToBeInterpolated      n x 2 matrix filled with x and y
                                  coordinates of the points the Lagrange
                                  polynomial must go through. The order of
                                  the polynomial depends on n.
                                  Currently, n can be either 3 or 4.

    Output:
      lagrangeCoefficient        n x 1 vector filled with the coefficients of
                                 the Lagrange polynomial. The first
                                 coefficient corresponds to the highest order
                                 term of the Lagrange polynomial.

    Changelog
    Date        Name                Comments
    29/09/2016  M. Van den Broeck   Creation.
    09/12/2016  M. Van den Broeck   Added third order Lagrangian interpolation

    Michael Van den Broeck ----- 25/08/2016
*/

    /// Check input
    // Currently, only 3 points are allowed
    Eigen::Vector2i allowedSize_order2;
    allowedSize_order2 << 3, 2;

    Eigen::Vector2i allowedSize_order3;
    allowedSize_order3 << 4, 2;

    if( (pointsToBeInterpolated_.rows() == allowedSize_order2( 0 ))
            && (pointsToBeInterpolated_.cols() == allowedSize_order2( 1 )) )
    {
        // Define x and y coordinates
        Eigen::Vector3d x,y;
        x = pointsToBeInterpolated_.col(0);
        y = pointsToBeInterpolated_.col(1);

        // Simplify notation so that it matches mathematical description
        double x0, x1, x2;
        x0 = x(0);
        x1 = x(1);
        x2 = x(2);
        double y0, y1, y2;
        y0 = y(0);
        y1 = y(1);
        y2 = y(2);

        // Calculation terms t
        double t01, t02, t11, t12, t21, t22;
        t01 = - (x1 + x2);
        t02 = x1*x2;

        t11 = - (x0 + x2);
        t12 = x0*x2;

        t21 = - (x0 + x1);
        t22 = x0*x1;

        // Calculate denominators
        double d0, d1, d2;
        d0 = (x0-x1) * (x0-x2);
        d1 = (x1-x0) * (x1-x2);
        d2 = (x2-x0) * (x2-x1);

        // Finally calculate L(x)
        double T1, T2, T3;
        T1 = y0/d0 + y1/d1 + y2/d2;
        T2 = y0*t01/d0 + y1*t11/d1 + y2*t21/d2;
        T3 = y0*t02/d0 + y1*t12/d1 + y2*t22/d2;

        // Give output
        Eigen::Matrix< double, 3, 1> lagrangeCoefficients;
        lagrangeCoefficients << T1,
                                T2,
                                T3;
        return lagrangeCoefficients;

    }
    else if( (pointsToBeInterpolated_.rows() == allowedSize_order3( 0 ))
            && (pointsToBeInterpolated_.cols() == allowedSize_order3( 1 )) )
    {
        // Define x and y coordinates
        Eigen::Vector4d x,y;
        x = pointsToBeInterpolated_.col(0);
        y = pointsToBeInterpolated_.col(1);

        // Simplify notation so that it matches mathematical description
        double x0, x1, x2, x3;
        x0 = x(0);
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        double y0, y1, y2, y3;
        y0 = y(0);
        y1 = y(1);
        y2 = y(2);
        y3 = y(3);

        // Calculation terms t
        double t01, t02, t03, t11, t12, t13, t21, t22, t23, t31, t32, t33;
        t01 = - (x1 + x2 + x3);
        t02 = x1*x2 + x1*x3 + x2*x3;
        t03 = - (x1*x2*x3);

        t11 = - (x0 + x2 + x3);
        t12 = x0*x2 + x0*x3 + x2*x3;
        t13 = - (x0*x2*x3);

        t21 = - (x0 + x1 + x3);
        t22 = x0*x1 + x0*x3 + x1*x3;
        t23 = - (x0*x1*x3);

        t31 = - (x0 + x1 + x2);
        t32 = x0*x1 + x0*x2 + x1*x2;
        t33 = - (x0*x1*x2);

        // Calculate denominators
        double d0, d1, d2, d3;
        d0 = (x0-x1) * (x0-x2) * (x0-x3);
        d1 = (x1-x0) * (x1-x2) * (x1-x3);
        d2 = (x2-x0) * (x2-x1) * (x2-x3);
        d3 = (x3-x0) * (x3-x1) * (x3-x2);

        // Finally calculate L(x)
        double T1, T2, T3, T4;
        T1 = y0/d0 + y1/d1 + y2/d2 + y3/d3;
        T2 = y0*t01/d0 + y1*t11/d1 + y2*t21/d2 + y3*t31/d3;
        T3 = y0*t02/d0 + y1*t12/d1 + y2*t22/d2 + y3*t32/d3;
        T4 = y0*t03/d0 + y1*t13/d1 + y2*t23/d2 + y3*t33/d3;

        // Give output
        Eigen::Matrix< double, 4, 1> lagrangeCoefficients;
        lagrangeCoefficients << T1,
                                T2,
                                T3,
                                T4;
        return lagrangeCoefficients;

    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Incorrect input size. Size of input should either be [3,2] or [4,2]." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

}



LagrangeInterpolator::~LagrangeInterpolator()
{

}

