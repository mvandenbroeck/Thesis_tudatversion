set(TEST   "${SRCROOT}/recurrenceRelations.h"
           "${SRCROOT}/problemRecurrenceRelations.cpp"
           "${SRCROOT}/spacecraft.cpp"
           "${SRCROOT}/stepSizeControlTSI.cpp"
           "${SRCROOT}/integrationSettings.cpp"
           "${SRCROOT}/lagrangeInterpolator.cpp"
           "${SRCROOT}/discreteForceModel.cpp"
           "${SRCROOT}/frameTransformation.cpp"
           "${SRCROOT}/taylorSeriesIntegrator.cpp"
)

message(STATUS "Source route is " ${SRCROOT})

# Add library.
add_library(unitTest_library STATIC ${TEST})
setup_library_target(unitTest_library "${CODEROOT}/lib")

# Add unitTestRecurrenceRelation application.
add_executable(application_unitTestRecurrenceRelations_tudatversion "${SRCROOT}/unitTests/unitTestProblemRecurrenceRelations.cpp")
setup_executable_target(application_unitTestRecurrenceRelations_tudatversion "${SRCROOT}")
target_link_libraries(application_unitTestRecurrenceRelations_tudatversion unitTest_library ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

# Add unitTestProblemRecurrenceRelation application.
add_executable(application_unitTestProblemRecurrenceRelations_tudatversion "${SRCROOT}/unitTests/unitTestProblemRecurrenceRelations.cpp")
setup_executable_target(application_unitTestProblemRecurrenceRelations_tudatversion "${SRCROOT}")
target_link_libraries(application_unitTestProblemRecurrenceRelations_tudatversion unitTest_library ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

# Add unitTestStepSizeControlTSI application.
add_executable(application_unitTestStepSizeControlTSI_tudatversion "${SRCROOT}/unitTests/unitTestStepSizeControlTSI.cpp")
setup_executable_target(application_unitTestStepSizeControlTSI_tudatversion "${SRCROOT}")
target_link_libraries(application_unitTestStepSizeControlTSI_tudatversion unitTest_library ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

# Add unitTestLagrangeInterpolator application.
add_executable(application_unitTestLagrangeInterpolator_tudatversion "${SRCROOT}/unitTests/unitTestLagrangeInterpolator.cpp")
setup_executable_target(application_unitTestLagrangeInterpolator_tudatversion "${SRCROOT}")
target_link_libraries(application_unitTestLagrangeInterpolator_tudatversion unitTest_library ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

# Add unitTestTaylorSeriesIntegrator application.
add_executable(application_unitTestTaylorSeriesIntegrator_tudatversion "${SRCROOT}/unitTests/unitTestTaylorSeriesIntegrator.cpp")
setup_executable_target(application_unitTestTaylorSeriesIntegrator_tudatversion "${SRCROOT}")
target_link_libraries(application_unitTestTaylorSeriesIntegrator_tudatversion unitTest_library tudat_input_output
tudat_basic_mathematics ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

# Add unitTestFrameTransformation application.
add_executable(application_unitTestFrameTransformation_tudatversion "${SRCROOT}/unitTests/unitTestFrameTransformation.cpp")
setup_executable_target(application_unitTestFrameTransformation_tudatversion "${SRCROOT}")
target_link_libraries(application_unitTestFrameTransformation_tudatversion unitTest_library tudat_input_output
tudat_basic_mathematics ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
