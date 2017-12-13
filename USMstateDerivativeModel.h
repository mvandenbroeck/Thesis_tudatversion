#ifndef USMSTATEDERIVATIVEMODEL_H
#define USMSTATEDERIVATIVEMODEL_H

#include "spacecraftGen.h"
#include "constants.h"

#include <boost/shared_ptr.hpp>


class USMstateDerivativeModel
{
public:
    USMstateDerivativeModel(
            SpacecraftGenPointer spacecraftGenPointer,
            ConstantsPointer constantsPointer );

    void computeCurrentAccelerations( const double &currentTime );

    Eigen::VectorXd computeUnifiedStateModelStateDerivative( const double currentTime, const Eigen::VectorXd& currentState );

    ~USMstateDerivativeModel(){}

private:
    SpacecraftGenPointer spacecraftGenPointer_;
    ConstantsPointer constantsPointer_;
};

//typedef boost::shared_ptr< USMstateDerivativeModel > USMstateDerivativeModelPointer;

#endif // USMSTATEDERIVATIVEMODEL_H
