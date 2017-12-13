#ifndef FRAMETRANSFORMATION_H
#define FRAMETRANSFORMATION_H
#include "Eigen/Core"

class FrameTransformation
{
public:
    FrameTransformation();

    Eigen::Matrix< double, 1, 3 > velocityFrameToUSMFrame( Eigen::Matrix< double, 1, 3> thrustInVelocityFrame,
                                                           double flightPathAngle );

    Eigen::Matrix< double, 1, 3 > velocityFrameToCartFrame(Eigen::Matrix< double, 1, 3 > thrustInVelocityFrame,
            Eigen::Matrix< double, 6, 1 > vehicleStateCartesian );

    ~FrameTransformation();
};

#endif // FRAMETRANSFORMATION_H
