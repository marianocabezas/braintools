#ifndef BRAINREGISTRATION_H
#define BRAINREGISTRATION_H

#include <vector>
#include <cstdlib>

#include <itkImageRegistrationMethod.h>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>

#include <itkVersorRigid3DTransform.h>
#include <itkAffineTransform.h>
#include <itkBSplineDeformableTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkScaleTransform.h>
#include <itkDemonsRegistrationFilter.h>


#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkVersorRigid3DTransformOptimizer.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>

#include <itkMutualInformationImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkMutualInformationHistogramImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkGradientDifferenceImageToImageMetric.h>
#include <itkKullbackLeiblerCompareHistogramImageToImageMetric.h>
#include <itkKappaStatisticImageToImageMetric.h>
#include <itkCorrelationCoefficientHistogramImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMeanSquaresHistogramImageToImageMetric.h>
#include <itkMeanReciprocalSquareDifferenceImageToImageMetric.h>

#include <itkBSplineResampleImageFunction.h>
#include <itkBSplineDecompositionImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkWarpImageFilter.h>

#include <itkCastImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkImageMomentsCalculator.h>
#include <itkSubtractImageFilter.h>
#include <itkDeformationFieldJacobianDeterminantFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkGradientImageFilter.h>

#include <imagedefinitions.h>
#include <registrationlevelupdate.h>
#include <brainio.h>
#include <brainsegmentation.h>


class BrainRegistration
{

    typedef double CoordinateRepType;

    // Registration methods
    typedef itk::ImageRegistrationMethod<
            SignedImageType,
            SignedImageType > RegistrationType;
	typedef itk::ImageRegistrationMethod<
            ProbabilityImageType,
            ProbabilityImageType > PRegistrationType;
    typedef itk::MultiResolutionImageRegistrationMethod<
            SignedImageType,
            SignedImageType > MultiResolutionRegistrationType;
    typedef itk::MultiResolutionImageRegistrationMethod<
            ProbabilityImageType,
            ProbabilityImageType > PMultiResolutionRegistrationType;
	typedef itk::DemonsRegistrationFilter<
			ProbabilityImageType,
			ProbabilityImageType,
			DeformationFieldType > DemonsRegistrationFilterType;
	typedef itk::MultiResolutionPDEDeformableRegistration<
			ProbabilityImageType, 
			ProbabilityImageType, 
			DeformationFieldType > MultiDemonsRegistrationFilterType;


    // Resamplers
    typedef itk::ResampleImageFilter<
            IntensityImageType,
            IntensityImageType > ResampleImageFilterType;
    typedef itk::ResampleImageFilter<
            SignedImageType,
            SignedImageType > ResampleSignedImageFilterType;
    typedef itk::ResampleImageFilter<
            ProbabilityImageType,
            ProbabilityImageType > ResampleAtlasFilterType;
    typedef itk::ResampleImageFilter<
            MaskImageType,
            MaskImageType > ResampleMaskFilterType;
    typedef itk::ResampleImageFilter<
            ResultsImageType,
            ResultsImageType > ResampleResultsFilterType;
    typedef itk::MultiResolutionPyramidImageFilter<
            SignedImageType,
            SignedImageType > MultiResolutionPyramidType;
    typedef itk::MultiResolutionPyramidImageFilter<
            ProbabilityImageType,
            ProbabilityImageType > PMultiResolutionPyramidType;
	typedef itk::WarpImageFilter<
			ProbabilityImageType,
			ProbabilityImageType,
			DeformationFieldType > WarperType;
	typedef itk::WarpImageFilter<
			MaskImageType,
			MaskImageType,
			DeformationFieldType > WarperMaskType;


    // Transform definitions
    typedef itk::MinimumMaximumImageCalculator< IntensityImageType > MinimumMaximumImageCalculatorType;
    typedef itk::MinimumMaximumImageCalculator< ProbabilityImageType > PMinimumMaximumImageCalculatorType;
    typedef itk::IdentityTransform< double, 3 > IdentityTransformType;
    typedef itk::ScaleTransform< double, 3 > ScaleTransformType;
    typedef itk::VersorRigid3DTransform< double > RigidTransformType;
	typedef RigidTransformType::VersorType VersorType;
	typedef VersorType::VectorType VectorType;
    typedef itk::AffineTransform< double, 3 > AffineTransformType;
    typedef itk::BSplineDeformableTransform< double, 3, 3 >     DeformableTransformType;
    typedef itk::CenteredTransformInitializer< RigidTransformType, IntensityImageType, IntensityImageType > TransformInitializerType;
    typedef itk::CenteredTransformInitializer< RigidTransformType, ProbabilityImageType, ProbabilityImageType > PTransformInitializerType;
    typedef DeformableTransformType::RegionType RegionType;
    typedef DeformableTransformType::SpacingType SpacingType;
    typedef DeformableTransformType::OriginType OriginType;
    typedef DeformableTransformType::ParametersType ParametersType;
    typedef DeformableTransformType::ImageType ParametersImageType;
    typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
    typedef itk::BSplineResampleImageFunction<ParametersImageType,CoordinateRepType> FunctionType;
    typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType> DecompositionType;
    typedef itk::ImageRegionIterator<ParametersImageType> Iterator;

    // Casting
    typedef itk::CastImageFilter< IntensityImageType, SignedImageType> IntensityToSignedCastType;

    // Optimizer
    typedef itk::VersorRigid3DTransformOptimizer VersorOptimizerType;
    typedef VersorOptimizerType::ScalesType VersorScalesType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef OptimizerType::ScalesType OptimizerScalesType;

    // Interpolator
    typedef itk::LinearInterpolateImageFunction<SignedImageType,double > InterpolatorType;
    typedef itk::LinearInterpolateImageFunction<ProbabilityImageType,double > PInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<MaskImageType,double> NNInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<ResultsImageType,double> NNAnatomicInterpolatorType;
    typedef itk::BSplineInterpolateImageFunction<IntensityImageType,double> BSplineImageInterpolatorType;
    typedef itk::BSplineInterpolateImageFunction<ProbabilityImageType,double> BSplineAtlasInterpolatorType;
    typedef itk::BSplineInterpolateImageFunction<SignedImageType,double> BSplineInterpolatorType;


    // Metric definitions
    typedef itk::ImageToImageMetric<SignedImageType,SignedImageType> GenericMetricType;
    typedef itk::MutualInformationImageToImageMetric<SignedImageType,SignedImageType> MIMetricType;
    typedef itk::MattesMutualInformationImageToImageMetric<SignedImageType,SignedImageType > MattesMIMetricType;
    typedef itk::MattesMutualInformationImageToImageMetric<ProbabilityImageType,ProbabilityImageType > PMattesMIMetricType;
    typedef itk::MutualInformationHistogramImageToImageMetric<SignedImageType,SignedImageType > MIHistogramMetricType;
    typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<SignedImageType,SignedImageType > NMIMetricType;
    typedef itk::GradientDifferenceImageToImageMetric<SignedImageType,SignedImageType> GDMetricType;
    typedef itk::KullbackLeiblerCompareHistogramImageToImageMetric<SignedImageType, SignedImageType> KLHistogramMetricType;
    typedef itk::KappaStatisticImageToImageMetric<SignedImageType, SignedImageType> KappaMetricType;
    typedef itk::CorrelationCoefficientHistogramImageToImageMetric<SignedImageType, SignedImageType> CorrelationMetricType;
    typedef itk::NormalizedCorrelationImageToImageMetric<SignedImageType, SignedImageType> NormalizedCorrelationMetricType;
    typedef itk::MeanSquaresImageToImageMetric<SignedImageType, SignedImageType> MSMetricType;
    typedef itk::MeanSquaresHistogramImageToImageMetric<SignedImageType, SignedImageType> MSHistogramMetricType;
    typedef itk::MeanReciprocalSquareDifferenceImageToImageMetric<SignedImageType, SignedImageType> MSReciprocalMetricType;

    // Observers
    typedef RegistrationLevelUpdate<MultiResolutionRegistrationType> ObserverType;

    // PreProcessing
    typedef itk::HistogramMatchingImageFilter< SignedImageType, SignedImageType > MatchingFilterType;
    typedef itk::RecursiveGaussianImageFilter< SignedImageType, SignedImageType > GaussianFilterType;

	// Others
	typedef itk::ImageMomentsCalculator < ProbabilityImageType >	BaseImageCalculatorType;
	typedef itk::SubtractImageFilter< ProbabilityImageType, ProbabilityImageType, ProbabilityImageType > SubtractionFilterType;
	typedef itk::DeformationFieldJacobianDeterminantFilter< DeformationFieldType >  JacobianFilterType;
	typedef itk::AddImageFilter< DeformationFieldType, DeformationFieldType > AddFilterType;
	typedef itk::MultiplyByConstantImageFilter< DeformationFieldType, ProbabilityPixelType, DeformationFieldType > MultiplyConstantFilterType;
	typedef itk::GradientImageFilter< ProbabilityImageType, float > GradientFilterType;
	typedef GradientFilterType::OutputImageType GradientImageType; 
	typedef GradientImageType::Pointer GradientImage;
	typedef itk::ImageRegionIterator<GradientImageType> GradientIterator;



public:
	/* Class methods */
    BrainRegistration();

    // Setters
    void SetFixed(IntensityImage fix);
    void SetMoving(IntensityImage mov);
    void SetMask(MaskImage mask);
    void SetAnatomic(ResultsImage mask);
    void AddAtlas(ProbabilityImage atlas);
    void SetAtlases(std::vector<ProbabilityImage> atlases);
    void ClearAtlases();
    // Getters
    IntensityImage GetFixed();
    IntensityImage GetMoving();
    MaskImage GetMask();
    ResultsImage GetAnatomic();
    std::vector<ProbabilityImage> GetAtlases();
    ProbabilityImage GetSimilarity(unsigned long radiusSize=1);

    // Registration methods
    void RigidRegistration(int levels = 3, int steps = 50);
    void AffineRegistration(int levels = 3, int steps = 50);
    void BSplineMultiRegistration(int levels = 2, unsigned int ngridnodesdim = 5, int steps = 200);
    void BSplineRegistration(unsigned int ngridnodesdim = 10, int steps = 200);

	/* Static methods */
	// registration methods
	static RigidTransform RigidRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int steps = 50);
	static AffineTransform AffineRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int steps = 50);
	static RigidTransform MultiRigidRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels = 2, int steps = 50);
    static AffineTransform MultiAffineRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels = 2, int steps = 50);
	static DeformableTransform BSplineMultiRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, AffineTransform affinetransf, int levels = 3, unsigned int ngridnodesdim = 8, int steps = 100);
	static DeformableTransform BSplineMultiRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels = 3, unsigned int ngridnodesdim = 8, int steps = 100);
	static DeformationField DemonsRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int steps = 50);
	static DeformationField MultiDemonsRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels = 2, int steps = 50);
	static ProbabilityImage Warp(ProbabilityImage fiximg, ProbabilityImage movimg, DeformationField deformation); 
	static MaskImage Warp(MaskImage fiximg, MaskImage movimg, DeformationField deformation);
	static DeformationField Sum(std::vector<DeformationField> deformations);
	static ProbabilityImage Jacobian(DeformationField deformation);
	static ProbabilityImage Divergence(DeformationField deformation);
	// Similarity
	template <class TTransform>
	static ProbabilityImage ResampleImage(ProbabilityImage fixed, ProbabilityImage moving, typename TTransform::Pointer transf) {
		/* Init */
		std::cout << "\tBrainRegistration::ResampleImage" << std::endl;
		std::cout << "\t\\- Creating image" << std::endl;
		ResampleAtlasFilterType::Pointer imresample;
		BSplineAtlasInterpolatorType::Pointer interpolator = BSplineAtlasInterpolatorType::New();

		// We resample the image to the rigid transform
		imresample = ResampleAtlasFilterType::New();
		imresample->SetInput(moving);
		imresample->SetSize(fixed->GetLargestPossibleRegion().GetSize());
		imresample->SetOutputOrigin(fixed->GetOrigin());
		imresample->SetOutputSpacing(fixed->GetSpacing());
		imresample->SetOutputDirection(fixed->GetDirection());
		imresample->SetInterpolator(interpolator);
		imresample->SetDefaultPixelValue( 0 );
		imresample->SetTransform(transf);
		try {
			imresample->Update();
		}
		catch( itk::ExceptionObject & err ) {
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}

		// We prepare the nighborhood iterators for both images
		return(imresample->GetOutput());
	}
	template <class TTransform>
	static IntensityImage ResampleImage(IntensityImage fixed, IntensityImage moving, typename TTransform::Pointer transf) {
		/* Init */
		std::cout << "\tBrainRegistration::ResampleImage" << std::endl;
		std::cout << "\t\\- Creating image" << std::endl;
		ResampleImageFilterType::Pointer imresample;
		BSplineImageInterpolatorType::Pointer interpolator = BSplineImageInterpolatorType::New();

		// We resample the image to the rigid transform
		imresample = ResampleImageFilterType::New();
		imresample->SetInput(moving);
		imresample->SetSize(fixed->GetLargestPossibleRegion().GetSize());
		imresample->SetOutputOrigin(fixed->GetOrigin());
		imresample->SetOutputSpacing(fixed->GetSpacing());
		imresample->SetOutputDirection(fixed->GetDirection());
		imresample->SetInterpolator(interpolator);
		imresample->SetDefaultPixelValue( 0 );
		imresample->SetTransform(transf);
		try {
			imresample->Update();
		}
		catch( itk::ExceptionObject & err ) {
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}

		// We prepare the neighborhood iterators for both images
		return(imresample->GetOutput());
	}
	template <class TTransform>
	static MaskImage ResampleImage(MaskImage fixed, MaskImage moving, typename TTransform::Pointer transf) {
		/* Init */
		std::cout << "\tBrainRegistration::ResampleImage" << std::endl;
		std::cout << "\t\\- Creating image" << std::endl;
		ResampleMaskFilterType::Pointer imresample;
		NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();

		// We resample the image to the rigid transform
		imresample = ResampleMaskFilterType::New();
		imresample->SetInput(moving);
		imresample->SetSize(fixed->GetLargestPossibleRegion().GetSize());
		imresample->SetOutputOrigin(fixed->GetOrigin());
		imresample->SetOutputSpacing(fixed->GetSpacing());
		imresample->SetOutputDirection(fixed->GetDirection());
		imresample->SetInterpolator(interpolator);
		imresample->SetTransform(transf);
		try {
			imresample->Update();
		}
		catch( itk::ExceptionObject & err ) {
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}

		// We prepare the nighborhood iterators for both images
		return(imresample->GetOutput());
	}
	template <class TTransform>
	static ResultsImage ResampleImage(ResultsImage fixed, ResultsImage moving, typename TTransform::Pointer transf) {
		/* Init */
		std::cout << "\tBrainRegistration::ResampleImage" << std::endl;
		std::cout << "\t\\- Creating image" << std::endl;
		ResampleResultsFilterType::Pointer imresample;
		NNAnatomicInterpolatorType::Pointer interpolator = NNAnatomicInterpolatorType::New();

		// We resample the image to the rigid transform
		imresample = ResampleResultsFilterType::New();
		imresample->SetInput(moving);
		imresample->SetSize(fixed->GetLargestPossibleRegion().GetSize());
		imresample->SetOutputOrigin(fixed->GetOrigin());
		imresample->SetOutputSpacing(fixed->GetSpacing());
		imresample->SetOutputDirection(fixed->GetDirection());
		imresample->SetInterpolator(interpolator);
		imresample->SetTransform(transf);
		try {
			imresample->Update();
		}
		catch( itk::ExceptionObject & err ) {
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}

		// We prepare the nighborhood iterators for both images
		return(imresample->GetOutput());
	}
	template <class TTransform>
	static ProbabilityImage GetSimilarity(ProbabilityImage fixed, ProbabilityImage moving, typename TTransform::Pointer affine, unsigned long radiusSize=1) {
		/* Init */
		BrainIO* brainio = new BrainIO();
		unsigned int i, j;
		std::cout << "\tBrainRegistration::GetSimilarity" << std::endl;
		std::cout << "\t\\- Creating image" << std::endl;
		ProbabilityImage similarity = ProbabilityImageType::New();
		similarity->SetRegions( fixed->GetLargestPossibleRegion() );
		similarity->SetSpacing( fixed->GetSpacing() );
		similarity->SetOrigin( fixed->GetOrigin() );
		similarity->SetDirection( fixed->GetDirection() );
		similarity->Allocate();
		similarity->FillBuffer(0);
		similarity->Update();

		// We resample the image to the last transform
		ProbabilityImage moved = ResampleImage<TTransform>(fixed, moving, affine);

		// We prepare the neighborhood iterators for both images
		std::cout << "\t\\- Preparing the iterators" << std::endl;
		ProbabilityNeighborhoodIterator::RadiusType radius;
		for (i = 0; i < ProbabilityImageType::ImageDimension; ++i)
			radius[i] = radiusSize;
		ProbabilityNeighborhoodIterator movIt = ProbabilityNeighborhoodIterator(radius,moved,moved->GetLargestPossibleRegion());
		ProbabilityNeighborhoodIterator fixIt = ProbabilityNeighborhoodIterator(radius,fixed,fixed->GetLargestPossibleRegion());

		std::cout << "\t\\- Cross-correlation computation" << std::endl;
		ProbabilityIterator simIt = ProbabilityIterator(similarity,similarity->GetLargestPossibleRegion());
		ProbabilityPixelType movMean, fixMean, ncc, fixStd, movStd;

		for (simIt.GoToBegin();!simIt.IsAtEnd();++simIt, ++movIt, ++fixIt) {
			if (movIt.GetCenterPixel()>0 && fixIt.GetCenterPixel()>0) {
				// Mean computation
				movMean = 0;
				fixMean = 0;
				movStd = 0;
				fixStd = 0;
				ncc = 0;
				for (j=0; j<movIt.Size(); j++) {
					movMean += movIt.GetPixel(j);
					fixMean += fixIt.GetPixel(j);
				}
				movMean = movMean/movIt.Size();
				fixMean = fixMean/fixIt.Size();

				// Normalised cross correlation computation
				for (j=0; j<movIt.Size(); j++) {
					fixStd += (fixIt.GetPixel(j) - fixMean)*(fixIt.GetPixel(j) - fixMean);
					movStd += (movIt.GetPixel(j) - movMean)*(movIt.GetPixel(j) - movMean);
					ncc += (movIt.GetPixel(j) - movMean)*(fixIt.GetPixel(j) - fixMean);
				}

				simIt.Set(fabs(ncc/(sqrt(fixStd)*sqrt(movStd))));
			}
			else
				simIt.Set(0);
		}

		return(similarity);
	}
	static ProbabilityImage Subtraction(ProbabilityImage minuend, ProbabilityImage subtrahend, RigidTransform transformation);

private:

    void ResampleBSpline();

    // Registration done flags
    bool rigid;
    bool affine;
    bool bspline;

    // Registration parameters
    int rigidlevels;
    int affinelevels;
    double bsplinenodes;

    // Image based pointers
    IntensityImage fiximg;
    IntensityRegionType fixreg;
    IntensitySizeType fixsize;
    OriginType origin;
    SpacingType spacing;
    IntensityDirectionType griddir;
    IntensityImage movimg;
    MaskImage maskimg;
    ResultsImage anaimg;
    IntensityRegionType movreg;
    IntensitySizeType movsize;
    OriginType movorigin;
    SpacingType movspacing;
    std::vector<ProbabilityImage> movatlas;

    // Region based pointers
    RegionType::SizeType gridsizeimg;
    RegionType::SizeType gridbordersize;
    RegionType::SizeType totalgridsize;

    // Transform pointers
    IdentityTransformType::Pointer idtransf;
    RigidTransformType::Pointer rigidtransf;
    AffineTransformType::Pointer affinetransf;
    DeformableTransformType::Pointer bsplinetransf;
    RegionType bsplinereg;
    TransformInitializerType::Pointer initializer;

    // Metric pointers
    MattesMIMetricType::Pointer mattesmi_metric;


    // Optimizer pointers
    OptimizerType::Pointer optimizer;
    VersorOptimizerType::Pointer versor_optimizer;

    // Interpolator pointers
    InterpolatorType::Pointer interpolator;
    BSplineImageInterpolatorType::Pointer bintinterpolator;
    BSplineAtlasInterpolatorType::Pointer bprinterpolator;
    BSplineInterpolatorType::Pointer binterpolator;
    NNInterpolatorType::Pointer ninterpolator;
    NNAnatomicInterpolatorType::Pointer nanainterpolator;

    // Registration methods
    RegistrationType::Pointer registration;
    MultiResolutionRegistrationType::Pointer multires_registration;

    // Resampling pointers
    ResampleImageFilterType::Pointer imresample;
    ResampleAtlasFilterType::Pointer atlasresample;
    ResampleMaskFilterType::Pointer maskresample;
    ResampleResultsFilterType::Pointer anatomicresample;
    MultiResolutionPyramidType::Pointer multires_pyramidfix;
    MultiResolutionPyramidType::Pointer multires_pyramidmov;

    // Casting pointers
    IntensityToSignedCastType::Pointer fixint2signed;
    IntensityToSignedCastType::Pointer movint2signed;
    IntensityToSignedCastType::Pointer transfint2signed;

    // Command pointer
    ObserverType::Pointer observer;

    // PreProcessing pointers
    MatchingFilterType::Pointer matcher;

    // Others
    MinimumMaximumImageCalculatorType::Pointer minmax;
};

#endif // BRAINREGISTRATION_H
