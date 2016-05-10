PROJECT(BrainToolsCARM)

SET(CMAKE_BUILD_TYPE "Release")
cmake_minimum_required(VERSION 2.5)

SET(ITK_DIR "/usr/local/ITK3/Build/")
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Please Set ITK_DIR")
ENDIF(ITK_FOUND)


SET (SEGMENTATION_SOURCES
    brainsegmentation.cpp
    gaussianestimator1d.cpp
    gaussianestimator2d.cpp
    gaussianestimator3d.cpp
    gaussianestimatornd.cpp
)
SET (SEGMENTATION_HEADERS
    brainsegmentation.h
    gaussianestimator.h
    gaussianestimator1d.h
    gaussianestimator2d.h
    gaussianestimator3d.h
    gaussianestimatornd.h
)
SET (SEGMENTATION_FILES
	${SEGMENTATION_SOURCES}
	${SEGMENTATION_HEADERS}
)

SET (REGISTRATION_SOURCES
    brainregistration.cpp
)
SET (REGISTRATION_HEADERS
    registrationlevelupdate.h
    brainregistration.h
)
SET (REGISTRATION_FILES
	${REGISTRATION_SOURCES}
	${REGISTRATION_HEADERS}
)

SET (PREPROCESSING_SOURCES
    itkN4MRIBiasFieldCorrectionImageFilter.txx
    itkBSplineControlPointImageFilter.txx
    brainpreprocessing.cpp
)
SET (PREPROCESSING_HEADERS
    itkN4MRIBiasFieldCorrectionImageFilter.h
    itkBSplineControlPointImageFilter.h
    brainpreprocessing.h
)
SET (PREPROCESSING_FILES
	${PREPROCESSING_SOURCES}
	${PREPROCESSING_HEADERS}
)

SET (OTHER_SOURCES
    brainio.cpp
)
SET (OTHER_HEADERS
    imagedefinitions.h
    brainio.h
)
SET (OTHER_FILES
	${OTHER_SOURCES}
	${OTHER_HEADERS}
)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)

SET(PRE_FILES ${PREPROCESSING_FILES} ${OTHER_FILES})
SET(COREG_FILES ${REGISTRATION_FILES} ${OTHER_FILES})
SET(ATLAS_FILES ${REGISTRATION_FILES} ${OTHER_FILES})
SET(TISSUE_FILES ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(SUB_FILES ${REGISTRATION_FILES} ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(DEFO_FILES ${PREPROCESSING_FILES} ${REGISTRATION_FILES} ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(POST_FILES ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(ANALYSE_FILES ${SEGMENTATION_FILES} ${OTHER_FILES})

ADD_EXECUTABLE( PreTool MACOSX_BUNDLE premain.cpp ${PRE_FILES})
ADD_EXECUTABLE( CoregTool MACOSX_BUNDLE coregmain.cpp ${COREG_FILES})
ADD_EXECUTABLE( AtlasTool MACOSX_BUNDLE atlasmain.cpp ${ATLAS_FILES})
ADD_EXECUTABLE( TissueTool MACOSX_BUNDLE tissuemain.cpp ${TISSUE_FILES})
ADD_EXECUTABLE( SubtractionTool MACOSX_BUNDLE submain.cpp ${SUB_FILES})
ADD_EXECUTABLE( DeformableTool MACOSX_BUNDLE defomain.cpp ${DEFO_FILES})
ADD_EXECUTABLE( PostTool MACOSX_BUNDLE postmain.cpp ${POST_FILES})
ADD_EXECUTABLE( AnalyseTool MACOSX_BUNDLE analysemain.cpp ${ANALYSE_FILES})
ADD_EXECUTABLE( HistogramTool MACOSX_BUNDLE histogrammain.cpp ${PRE_FILES})
ADD_EXECUTABLE( Test testmain.cpp )


TARGET_LINK_LIBRARIES( PreTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( CoregTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( AtlasTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( TissueTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( SubtractionTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( DeformableTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( PostTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( AnalyseTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)
TARGET_LINK_LIBRARIES( HistogramTool
        ITKCommon
        ITKIO
        ITKStatistics
        ITKAlgorithms
        ITKBasicFilters
)



