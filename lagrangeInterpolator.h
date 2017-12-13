#ifndef LAGRANGEINTERPOLATOR_H
#define LAGRANGEINTERPOLATOR_H

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

class LagrangeInterpolator
{
public:
    // Constructor
    LagrangeInterpolator( Eigen::Matrix< double, Eigen::Dynamic, 2 > pointsToBeInterpolated);


    Eigen::Matrix<double, Eigen::Dynamic, 1> computeCoefficientsOfInterpolatingPolynomial();

    ~LagrangeInterpolator();

private:
    Eigen::Matrix< double, Eigen::Dynamic, 2 > pointsToBeInterpolated_;
};

typedef boost::shared_ptr< LagrangeInterpolator > LagrangeInterpolatorPointer;

#endif // LAGRANGEINTERPOLATOR_H
