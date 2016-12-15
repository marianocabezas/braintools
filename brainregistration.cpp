#include <brainregistration.h>


BrainRegistration::BrainRegistration() {
    rigid = false;
    affine = false;
    bspline = false;

    rigidlevels = 0;
    affinelevels = 0;
    bsplinenodes = 0;

    // Transforms init
    idtransf = IdentityTransformType::New();
    rigidtransf = RigidTransformType::New();
    affinetransf = AffineTransformType::New();
    bsplinetransf = DeformableTransformType::New();
    initializer = TransformInitializerType::New();

    // Optimizer init
    versor_optimizer = VersorOptimizerType::New();
    optimizer = OptimizerType::New();

    // Interpolator init
    interpolator = InterpolatorType::New();
    binterpolator = BSplineInterpolatorType::New();
    bintinterpolator = BSplineImageInterpolatorType::New();
    bprinterpolator = BSplineAtlasInterpolatorType::New();
    ninterpolator = NNInterpolatorType::New();
    nanainterpolator = NNAnatomicInterpolatorType::New();

    // Metric init
    mattesmi_metric = MattesMIMetricType::New();
    mattesmi_metric->ComputeGradientOff();
    mattesmi_metric->SetUseCachingOfBSplineWeights(false);  

    // Registration method init
    registration = RegistrationType::New();
    registration->SetMetric(mattesmi_metric);
    registration->SetOptimizer(optimizer);
    registration->SetInterpolator(interpolator);
    multires_registration = MultiResolutionRegistrationType::New();
    multires_registration->SetMetric(mattesmi_metric);
    multires_registration->SetOptimizer(optimizer);
    multires_registration->SetInterpolator(interpolator);

    // Resampling init
    imresample = ResampleImageFilterType::New();
    imresample->SetTransform(this->idtransf);
    imresample->ReleaseDataFlagOn();
    atlasresample = ResampleAtlasFilterType::New();
    atlasresample->SetTransform(this->idtransf);
    atlasresample->ReleaseDataFlagOn();
    maskresample = ResampleMaskFilterType::New();
    maskresample->SetTransform(this->idtransf);
    maskresample->SetInterpolator(ninterpolator);
    maskresample->ReleaseDataFlagOn();
    anatomicresample = ResampleResultsFilterType::New();
    anatomicresample->SetTransform(this->idtransf);
    anatomicresample->SetInterpolator(nanainterpolator);
    anatomicresample->ReleaseDataFlagOn();
    multires_pyramidfix = MultiResolutionPyramidType::New();
    multires_registration->SetFixedImagePyramid(this->multires_pyramidfix);
    multires_pyramidmov = MultiResolutionPyramidType::New();
    multires_registration->SetMovingImagePyramid(this->multires_pyramidmov);

    // Observer init
    observer = ObserverType::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );
    versor_optimizer->AddObserver( itk::IterationEvent(), observer );
    multires_registration->AddObserver( itk::IterationEvent(), observer );

    // Casting init
    fixint2signed = IntensityToSignedCastType::New();
    movint2signed = IntensityToSignedCastType::New();
    transfint2signed = IntensityToSignedCastType::New();

    // Preprocessing init
    matcher = MatchingFilterType::New();

    //Others
    minmax = MinimumMaximumImageCalculatorType::New();
}

void BrainRegistration::SetFixed(IntensityImage fix) {

    // Setting fixed image
    this->fiximg=fix;

    // Casting to a signed type for registration
    fixint2signed->SetInput(this->fiximg);
    fixint2signed->Update();

    // Registration parameters
    fixreg = fiximg->GetBufferedRegion();
    fixsize = this->fixreg.GetSize();
    origin = fiximg->GetOrigin();
    spacing = fiximg->GetSpacing();
    griddir = fiximg->GetDirection();

    multires_registration->SetFixedImage(fixint2signed->GetOutput());
    registration->SetFixedImage(fixint2signed->GetOutput());
    matcher->SetReferenceImage(fixint2signed->GetOutput());
    initializer->SetFixedImage(fiximg);

    // Intensity thresholds
    minmax->SetImage(fiximg);
    minmax->Compute();
    //mattesmi_metric->SetFixedImageSamplesIntensityThreshold(minmax->GetMaximum()*0.05);


    // No registration done yet
    affine = false;
    rigid = false;
}

void BrainRegistration::SetMoving(IntensityImage mov) {
    // Setting moving image
    this->movimg=mov;

    // Retrieving image information
    movreg = this->fiximg->GetBufferedRegion();
    movsize = this->fixreg.GetSize();
    movorigin = this->fiximg->GetOrigin();
    movspacing = this->fiximg->GetSpacing();

    // Casting to a signed type for registration
    movint2signed->SetInput(this->movimg);
    movint2signed->Update();

    //SignedImageType::Pointer down_movimg = this->DownsampleImage(movint2signed->GetOutput());
    //multires_registration->SetMovingImage(down_movimg);
    //registration->SetMovingImage(down_movimg);
    multires_registration->SetMovingImage(movint2signed->GetOutput());
    registration->SetMovingImage(movint2signed->GetOutput());
    initializer->SetMovingImage(movimg);

    // No registration done yet
    affine = false;
    rigid = false;
}

void BrainRegistration::SetMask(MaskImage mask) {
    this->maskimg=mask;
    this->maskresample->SetInput(this->maskimg);
}

void BrainRegistration::SetAnatomic(ResultsImage anatomic) {
    this->anaimg=anatomic;
    this->anatomicresample->SetInput(this->anaimg);
}

void BrainRegistration::AddAtlas(ProbabilityImage atlas) {
    this->movatlas.push_back(atlas);
}
void BrainRegistration::SetAtlases(std::vector<ProbabilityImage> atlases) {
    this->movatlas.clear();
    for (unsigned int i=0;i<atlases.size();i++) {
        this->movatlas.push_back(atlases[i]);
    }
}

void BrainRegistration::ClearAtlases() {
    this->movatlas.clear();
}

void BrainRegistration::RigidRegistration(int levels, int steps) {
	std::cout << "\tBrainRegistration::RigidRegistration" << std::endl;
    if (!rigid || rigidlevels < levels) {
        this->bspline=false;
        this->affine=false;

        rigidlevels = levels;

        std::cout << "\t\\- Initializing rigid transform" << std::endl;
        // Registration initialization
		VersorType rotation;
		VectorType axis;
		axis[0] = 0.0;
		axis[1] = 0.0;
		axis[2] = 1.0;
		rotation.Set( axis, 0.0 );
		rigidtransf->SetRotation( rotation );
        initializer->SetTransform(rigidtransf);
		initializer->MomentsOn();
		//initializer->GeometryOn();
        initializer->InitializeTransform();

        // Optimizer parameters
        const double translscale = 1.0 / 1000.0;
        VersorScalesType optimizerscales(rigidtransf->GetNumberOfParameters());
        optimizerscales[0] = 1.0;
        optimizerscales[1] = 1.0;
        optimizerscales[2] = 1.0;
        optimizerscales[3] = translscale;
        optimizerscales[4] = translscale;
        optimizerscales[5] = translscale;
        versor_optimizer->SetScales( optimizerscales );
        versor_optimizer->SetMaximumStepLength( 0.2000  );
        versor_optimizer->SetMinimumStepLength( 0.0001 );
        versor_optimizer->SetNumberOfIterations( steps );
		mattesmi_metric->SetNumberOfHistogramBins(50);
		mattesmi_metric->UseAllPixelsOn();

        std::cout << "\t\\- Starting rigid transform" << std::endl;
        // Registration
        multires_registration->SetOptimizer(versor_optimizer);
        multires_registration->SetTransform(rigidtransf);
        multires_registration->SetFixedImageRegion(fixreg);
        multires_registration->SetNumberOfLevels(levels);
        multires_registration->SetInitialTransformParameters(rigidtransf->GetParameters());
        try {
            //registration->StartRegistration();
            multires_registration->StartRegistration();
            this->rigid = true;
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return;
        }
        rigidtransf->SetParameters(multires_registration->GetLastTransformParameters());
        std::cout << "\t\\- Rigid transform completed" << std::endl;
        std::cout << std::endl;
    }
    imresample->SetTransform(this->rigidtransf);
    atlasresample->SetTransform(this->rigidtransf);
    maskresample->SetTransform(this->rigidtransf);
    anatomicresample->SetTransform(this->rigidtransf);
}

void BrainRegistration::AffineRegistration(int levels, int steps) {
	std::cout << "\tBrainRegistration::AffineRegistration" << std::endl;
    if (!affine || affinelevels < levels) {
        this->bspline=false;

        affinelevels = levels;

        // Check if images are already rigidly registered
        if (!rigid) {
            this->RigidRegistration();
        }

        std::cout << "\t\\- Starting affine transform" << std::endl;
        // Transform parameters
        affinetransf->SetCenter(rigidtransf->GetCenter());
        affinetransf->SetTranslation(rigidtransf->GetTranslation());
        affinetransf->SetMatrix(rigidtransf->GetMatrix());

        // Optimizer parameters
        const double translscale = 1.0 / 1000.0;
        OptimizerScalesType optimizerscales(affinetransf->GetNumberOfParameters());
        optimizerscales[0] = 1.0;
        optimizerscales[1] = 1.0;
        optimizerscales[2] = 1.0;
        optimizerscales[3] = 1.0;
        optimizerscales[4] = 1.0;
        optimizerscales[5] = 1.0;
        optimizerscales[6] = 1.0;
        optimizerscales[7] = 1.0;
        optimizerscales[8] = 1.0;
        optimizerscales[9]  = translscale;
        optimizerscales[10] = translscale;
        optimizerscales[11] = translscale;
        optimizer->SetScales( optimizerscales );
        optimizer->SetMaximumStepLength( 0.2000  );
        optimizer->SetMinimumStepLength( 0.0001 );
        optimizer->SetNumberOfIterations( steps );

        // Metric parameters
        mattesmi_metric->SetNumberOfHistogramBins(50);
		mattesmi_metric->UseAllPixelsOn();

        // Registration
        multires_registration->SetOptimizer(optimizer);
        multires_registration->SetTransform(affinetransf);
        multires_registration->SetNumberOfLevels(levels);
        multires_registration->SetInitialTransformParameters(affinetransf->GetParameters());
        try {
            //registration->StartRegistration();
            multires_registration->StartRegistration();
            affine = true;
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return;
        }
        affinetransf->SetParameters(multires_registration->GetLastTransformParameters());
        std::cout << "\t\\- Affine transform completed" << std::endl;
        std::cout << std::endl;
    }
    imresample->SetTransform(this->affinetransf);
    atlasresample->SetTransform(this->affinetransf);
    maskresample->SetTransform(this->affinetransf);
    anatomicresample->SetTransform(this->affinetransf);
}

void BrainRegistration::BSplineMultiRegistration(int levels, unsigned int ngridnodesdim,int steps) {
	std::cout << "\tBrainRegistration::BSplineMultiRegistration" << std::endl;
    if (!bspline || this->bsplinenodes != ngridnodesdim) {
        this->bsplinenodes = ngridnodesdim;

        // Check if images are already affinely registered
        if (!affine) {
            this->AffineRegistration();
        }

        std::cout << "\t\\- Starting bspline transform" << std::endl;


        // Grid estimation
        gridsizeimg.Fill(ngridnodesdim);
        gridbordersize.Fill(3);    // Border for spline order = 3 ( 1 lower, 2 upper )
        totalgridsize = gridsizeimg + gridbordersize;
        this->bsplinereg.SetSize(totalgridsize);

        for(unsigned int r=0; r<3; r++) {
            spacing[r] *= static_cast<double>(fixsize[r] - 1)  /
                        static_cast<double>(gridsizeimg[r] - 1);
        }

        OriginType gridorig = origin - griddir * spacing;

        // Transform parameters
        bsplinetransf->SetGridSpacing(spacing);
        bsplinetransf->SetGridOrigin(gridorig);
        bsplinetransf->SetGridRegion(this->bsplinereg);
        bsplinetransf->SetGridDirection(griddir);
        bsplinetransf->SetBulkTransform(affinetransf);
        unsigned int nbsplineparam = bsplinetransf->GetNumberOfParameters();
        ParametersType initdeftransfparam(nbsplineparam);
        initdeftransfparam.Fill( 0.0 );
        bsplinetransf->SetParameters(initdeftransfparam);

        // Optimizer parameters
        OptimizerScalesType optimizerscales = OptimizerScalesType(nbsplineparam);
        optimizerscales.Fill( 1.0 );
        optimizer->SetScales( optimizerscales );
        optimizer->SetMaximumStepLength( 10.0 );
        optimizer->SetMinimumStepLength(  0.01 );
        optimizer->SetNumberOfIterations( steps );

        // Metric parameters
        mattesmi_metric->SetUseExplicitPDFDerivatives(false);
        try {
            mattesmi_metric->SetNumberOfSpatialSamples(static_cast<int>(fixreg.GetNumberOfPixels()*0.5));
        }
        catch( itk::ExceptionObject ) {}

        // Registration
        multires_registration->SetOptimizer(optimizer);
        multires_registration->SetInterpolator(interpolator);
        multires_registration->SetInitialTransformParameters(bsplinetransf->GetParameters() );
        multires_registration->SetNumberOfLevels(levels);
        multires_registration->SetTransform(bsplinetransf);
        try {
          multires_registration->StartRegistration();
          bspline = true;
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return;
        }
        bsplinetransf->SetParameters(multires_registration->GetLastTransformParameters());
        std::cout << "\t\\- Bspline transform completed" << std::endl;
        std::cout << std::endl;
    }

    //this->ResampleBSpline();

    imresample->SetInterpolator(bintinterpolator);
    imresample->SetTransform(this->bsplinetransf);
    atlasresample->SetInterpolator(bprinterpolator);
    atlasresample->SetTransform(this->bsplinetransf);
    maskresample->SetTransform(this->bsplinetransf);
    anatomicresample->SetTransform(this->bsplinetransf);
}

void BrainRegistration::BSplineRegistration(unsigned int ngridnodesdim,int steps) {
	std::cout << "\tBrainRegistration::BSplineRegistration" << std::endl;
    if (!bspline || this->bsplinenodes != ngridnodesdim) {
        this->bsplinenodes = ngridnodesdim;

        // Check if images are already affinely registered
        if (!affine) {
            this->AffineRegistration();
        }

        std::cout << std::endl;
        std::cout << "\t\\- Starting bspline transform" << std::endl;


        // Grid estimation
        gridsizeimg.Fill(ngridnodesdim);
        gridbordersize.Fill(3);    // Border for spline order = 3 ( 1 lower, 2 upper )
        totalgridsize = gridsizeimg + gridbordersize;
        this->bsplinereg.SetSize(totalgridsize);

        for(unsigned int r=0; r<3; r++) {
            spacing[r] *= static_cast<double>(fixsize[r] - 1)  /
                        static_cast<double>(gridsizeimg[r] - 1);
        }

        OriginType gridorig = origin - griddir * spacing;

        // Transform parameters
        bsplinetransf->SetGridSpacing(spacing);
        bsplinetransf->SetGridOrigin(gridorig);
        bsplinetransf->SetGridRegion(this->bsplinereg);
        bsplinetransf->SetGridDirection(griddir);
        bsplinetransf->SetBulkTransform(affinetransf);
        unsigned int nbsplineparam = bsplinetransf->GetNumberOfParameters();
        ParametersType initdeftransfparam(nbsplineparam);
        initdeftransfparam.Fill( 0.0 );
        bsplinetransf->SetParameters(initdeftransfparam);

        // Optimizer parameters
        OptimizerScalesType optimizerscales = OptimizerScalesType(nbsplineparam);
        optimizerscales.Fill( 1.0 );
        optimizer->SetScales( optimizerscales );
        optimizer->SetMaximumStepLength( 10.0 );
        optimizer->SetMinimumStepLength(  0.01 );
        optimizer->SetNumberOfIterations( steps );

        // Metric parameters
        mattesmi_metric->SetNumberOfSpatialSamples(fixreg.GetNumberOfPixels());

        // Registration
        registration->SetOptimizer(optimizer);
        registration->SetInterpolator(interpolator);
        registration->SetInitialTransformParameters(bsplinetransf->GetParameters() );
        registration->SetTransform(bsplinetransf);
        try {
          registration->StartRegistration();
          bspline = true;
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return;
        }
        bsplinetransf->SetParameters(registration->GetLastTransformParameters());
        std::cout << "\t\\- Bspline transform completed" << std::endl;
        std::cout << std::endl;
    }


    imresample->SetInterpolator(bintinterpolator);
    imresample->SetTransform(this->bsplinetransf);
    atlasresample->SetInterpolator(bprinterpolator);
    atlasresample->SetTransform(this->bsplinetransf);
    maskresample->SetTransform(this->bsplinetransf);
    anatomicresample->SetTransform(this->bsplinetransf);
}

IntensityImage BrainRegistration::GetMoving() {

    // We resample the image to the last transform
    imresample->SetInput(this->movimg);
    imresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    imresample->SetOutputOrigin(this->fiximg->GetOrigin());
    imresample->SetOutputSpacing(this->fiximg->GetSpacing());
    imresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        imresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return imresample->GetOutput();

}

MaskImage BrainRegistration::GetMask() {

    // We resample the image to the last transform
    maskresample->SetInput(this->maskimg);
    maskresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    maskresample->SetOutputOrigin(this->fiximg->GetOrigin());
    maskresample->SetOutputSpacing(this->fiximg->GetSpacing());
    maskresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        maskresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return maskresample->GetOutput();
}

ResultsImage BrainRegistration::GetAnatomic() {

    // We resample the image to the last transform
    anatomicresample->SetInput(this->anaimg);
    anatomicresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    anatomicresample->SetOutputOrigin(this->fiximg->GetOrigin());
    anatomicresample->SetOutputSpacing(this->fiximg->GetSpacing());
    anatomicresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        anatomicresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return anatomicresample->GetOutput();

}

IntensityImage BrainRegistration::GetFixed() {
    return(this->fiximg);
}

std::vector<ProbabilityImage> BrainRegistration::GetAtlases() {

    std::vector<ProbabilityImage> atlases;

    // We resample all probabilistic atlases to the last transform
    for (unsigned int i=0;i<this->movatlas.size();i++) {
        atlasresample = ResampleAtlasFilterType::New();
        if (bspline) atlasresample->SetTransform(this->bsplinetransf);
        else if (affine) atlasresample->SetTransform(this->affinetransf);
        else if (rigid) atlasresample->SetTransform(this->rigidtransf);
        else atlasresample->SetTransform(this->idtransf);
        atlasresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
        atlasresample->SetOutputOrigin(this->fiximg->GetOrigin());
        atlasresample->SetOutputSpacing(this->fiximg->GetSpacing());
        atlasresample->SetOutputDirection(this->fiximg->GetDirection());
        atlasresample->SetInput(this->movatlas[i]);
        try {
            atlasresample->Update();
            atlases.push_back(atlasresample->GetOutput());
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return atlases;
        }
    }

    return atlases;
}

void BrainRegistration::ResampleBSpline() {
    // Resampling B-spline transform
    unsigned int num_parameters,gridnodesxdir,dim,counter = 0;

    // We double the nodes per grid at each level as the resolution increases and change the B-Splines parameters
    RegionType gridreg;
    RegionType::SizeType gridsize = this->bsplinetransf->GetGridRegion().GetSize();
    gridnodesxdir = gridsize[0]-this->bsplinetransf->SplineOrder;
    gridsize.Fill(gridnodesxdir*2+this->bsplinetransf->SplineOrder);
    gridreg.SetSize(gridsize);

    // We also update the origin and the spacing of the new grid
    SpacingType spacing = fixint2signed->GetOutput()->GetSpacing();
    OriginType  origin  = fixint2signed->GetOutput()->GetOrigin();
    SignedDirectionType griddir = fixint2signed->GetOutput()->GetDirection();
    SignedSizeType fixedsize = fixint2signed->GetOutput()->GetBufferedRegion().GetSize();
    for(dim=0; dim<3; dim++) {
        spacing[dim] *= static_cast<double>(fixedsize[dim] - 1) / static_cast<double>(gridsize[dim] - 1);
    }
    origin -= griddir * spacing;

    // We finally set the number of parameters
    num_parameters = gridsize[0] * gridsize[1] * gridsize[2] * 3;
    ParametersType new_parameters( num_parameters );
    new_parameters.Fill( 0.0 );

    // We create the resampling objects
    ResamplerType::Pointer upsampler = ResamplerType::New();
    FunctionType::Pointer function = FunctionType::New();
    DecompositionType::Pointer decomposition = DecompositionType::New();

    ScaleTransformType::Pointer scale = ScaleTransformType::New();
    ScaleTransformType::ScaleType scale_vector;
    scale_vector.Fill(2.0);
    scale->SetScale(scale_vector);

    for ( dim = 0; dim < 3; dim++ ) {

        // We prepare the grid upsampling
        upsampler->SetInput( this->bsplinetransf->GetCoefficientImage()[dim] );
        upsampler->SetInterpolator( function );
        upsampler->SetTransform( scale );
        upsampler->SetSize( gridsize );
        upsampler->SetOutputSpacing( spacing );
        upsampler->SetOutputOrigin( origin );

        // We prepare the new grid
        decomposition->SetSplineOrder( this->bsplinetransf->SplineOrder );
        decomposition->SetInput( upsampler->GetOutput() );
        decomposition->Update();

        // We copy the new coefficients
        ParametersImageType::Pointer new_coefficients = decomposition->GetOutput();
        Iterator it( new_coefficients, this->bsplinetransf->GetGridRegion() );
        while ( !it.IsAtEnd() ) {
            new_parameters[ counter++ ] = it.Get()*2;
            ++it;
        }
    }

    // We set the final grid parameters
    this->bsplinetransf->SetGridRegion( gridreg );
    this->bsplinetransf->SetGridSpacing( spacing );
    this->bsplinetransf->SetGridOrigin( origin );
    this->bsplinetransf->SetParameters( new_parameters );

}

ProbabilityImage BrainRegistration::GetSimilarity(unsigned long radiusSize) {
    unsigned int i;
    ProbabilityImage similarity = ProbabilityImageType::New();
    similarity->SetRegions( this->fiximg->GetLargestPossibleRegion() );
    similarity->SetSpacing( fiximg->GetSpacing() );
    similarity->SetOrigin( fiximg->GetOrigin() );
    similarity->SetDirection( fiximg->GetDirection() );
    similarity->Allocate();
    similarity->FillBuffer(0);
    similarity->Update();

    // We resample the image to the last transform
    imresample->SetInput(this->movimg);
    imresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    imresample->SetOutputOrigin(this->fiximg->GetOrigin());
    imresample->SetOutputSpacing(this->fiximg->GetSpacing());
    imresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        imresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    // We prepare the nighborhood iterators for both images
    IntensityImage movedimg = imresample->GetOutput();
    IntensityNeighborhoodIterator::RadiusType radius;
    for (i = 0; i < IntensityImageType::ImageDimension; ++i)
        radius[i] = radiusSize;
    IntensityNeighborhoodIterator movIt = IntensityNeighborhoodIterator(radius,movedimg,movedimg->GetRequestedRegion());
    IntensityNeighborhoodIterator fixIt = IntensityNeighborhoodIterator(radius,fiximg,fiximg->GetRequestedRegion());
    ProbabilityIterator simIt = ProbabilityIterator(similarity,similarity->GetRequestedRegion());
    ProbabilityPixelType movMean, fixMean, ncc, fixStd, movStd;

    for (movIt.GoToBegin();!movIt.IsAtEnd();++movIt,++fixIt,++simIt) {
        if (movIt.GetCenterPixel()>0 && fixIt.GetCenterPixel()>0) {
            // Mean computation
            movMean = 0;
            fixMean = 0;
            movStd = 0;
            fixStd = 0;
            ncc = 0;
            for (i=0; i<movIt.Size(); i++) {
                movMean += movIt.GetPixel(i);
                fixMean += fixIt.GetPixel(i);
            }
            movMean = movMean/movIt.Size();
            fixMean = fixMean/fixIt.Size();

            // Normalised cross correlation computation
            for (i=0; i<movIt.Size(); i++) {
                fixStd += (fixIt.GetPixel(i) - fixMean)*(fixIt.GetPixel(i) - fixMean);
                movStd += (movIt.GetPixel(i) - movMean)*(movIt.GetPixel(i) - movMean);
                ncc += (movIt.GetPixel(i) - movMean)*(fixIt.GetPixel(i) - fixMean);
            }

            simIt.Set(fabs(ncc/(sqrt(fixStd)*sqrt(movStd))));
        }
        else
            simIt.Set(0);
    }

    return(similarity);
}

RigidTransform BrainRegistration::RigidRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int steps) {
	std::cout << "\tBrainRegistration::RigidRegistration" << std::endl;
    /* Init */
    // Definitions
    const double translscale = 1.0 / 1000.0;
    ProbabilityRegionType fixreg = movimg->GetLargestPossibleRegion();
    // Observer
    ObserverType::Pointer observer = ObserverType::New();
    // Metric parameters
    PMattesMIMetricType::Pointer mattesmi_metric = PMattesMIMetricType::New();
	mattesmi_metric->SetNumberOfHistogramBins(50);
	mattesmi_metric->UseAllPixelsOn();
    // Transformations
	VersorType rotation;
	VectorType axis;
    RigidTransformType::Pointer rigidtransf = RigidTransformType::New();
    // Optimizers
    VersorOptimizerType::Pointer versor_optimizer = VersorOptimizerType::New();
    versor_optimizer->AddObserver( itk::IterationEvent(), observer );
    // Interpolator
    PInterpolatorType::Pointer interpolator = PInterpolatorType::New();
    // Registration framework
    PRegistrationType::Pointer registration = PRegistrationType::New();
    registration->SetMovingImage(movimg);
    registration->SetFixedImage(fiximg);
    registration->SetMetric(mattesmi_metric);
    registration->SetInterpolator(interpolator);
    registration->AddObserver( itk::IterationEvent(), observer );

    /* Rigid */
    std::cout << std::endl;
    std::cout << "\t\\-- Initializing rigid transform" << std::endl;
    // Registration initialization
	axis[0] = 0.0;
	axis[1] = 0.0;
	axis[2] = 1.0;
	rotation.Set( axis, 0.0 );
	rigidtransf->SetRotation( rotation );
    PTransformInitializerType::Pointer initializer = PTransformInitializerType::New();
    initializer->SetFixedImage(fiximg);
    initializer->SetMovingImage(movimg);
    initializer->SetTransform(rigidtransf);
    initializer->MomentsOn();
	//initializer->GeometryOn();
    initializer->InitializeTransform();

    // Optimizer parameters
    VersorScalesType versoroptimizerscales(rigidtransf->GetNumberOfParameters());
    versoroptimizerscales[0] = 1.0;
    versoroptimizerscales[1] = 1.0;
    versoroptimizerscales[2] = 1.0;
    versoroptimizerscales[3] = translscale;
    versoroptimizerscales[4] = translscale;
    versoroptimizerscales[5] = translscale;
    versor_optimizer->SetScales( versoroptimizerscales );
    versor_optimizer->SetMaximumStepLength( 2  );
    versor_optimizer->SetMinimumStepLength( 0.001 );
	versor_optimizer->SetRelaxationFactor(0.5);
    versor_optimizer->SetNumberOfIterations( steps );
	versor_optimizer->MinimizeOn();

    // Registration
    registration->SetOptimizer(versor_optimizer);
    registration->SetTransform(rigidtransf);
    registration->SetFixedImageRegion(fixreg);
    registration->SetInitialTransformParameters(rigidtransf->GetParameters());
    std::cout << "\t\\-- Starting rigid transform" << std::endl;
    try {
        registration->StartRegistration();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(rigidtransf);
    }
    rigidtransf->SetParameters(registration->GetLastTransformParameters());
    std::cout << "\t\\-- Rigid transform completed" << std::endl;
    std::cout << std::endl;

	return(rigidtransf);
}

AffineTransform BrainRegistration::AffineRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int steps) {
	std::cout << "\tBrainRegistration::AffineRegistration" << std::endl;
    /* Init */
    // Definitions
    const double translscale = 1.0 / 1000.0;
    IntensityRegionType fixreg = fiximg->GetLargestPossibleRegion();
    // Observer
    ObserverType::Pointer observer = ObserverType::New();
    // Metric parameters
    PMattesMIMetricType::Pointer mattesmi_metric = PMattesMIMetricType::New();
	mattesmi_metric->SetNumberOfHistogramBins(50);
	mattesmi_metric->UseAllPixelsOn();
    // Transformations
	RigidTransformType::Pointer rigidtransf;
    AffineTransformType::Pointer affinetransf = AffineTransformType::New();
    // Optimizers
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );
    // Interpolator
    PInterpolatorType::Pointer interpolator = PInterpolatorType::New();
    // Registration framework
    PRegistrationType::Pointer registration = PRegistrationType::New();
    registration->SetMovingImage(movimg);
    registration->SetFixedImage(fiximg);
    registration->SetMetric(mattesmi_metric);
    registration->SetInterpolator(interpolator);
    registration->AddObserver( itk::IterationEvent(), observer );

    /* Rigid */
    rigidtransf = RigidRegistration(fiximg, movimg, steps);

    /* Affine */
    std::cout << std::endl;
    std::cout << "\t\\-- Starting affine transform" << std::endl;
	std::cout << std::endl;
    // Transform parameters
    affinetransf->SetCenter(rigidtransf->GetCenter());
    affinetransf->SetTranslation(rigidtransf->GetTranslation());
    affinetransf->SetMatrix(rigidtransf->GetMatrix());

    // Optimizer parameters
    OptimizerScalesType optimizerscales(affinetransf->GetNumberOfParameters());
    optimizerscales[0] = 1.0;
    optimizerscales[1] = 1.0;
    optimizerscales[2] = 1.0;
    optimizerscales[3] = 1.0;
    optimizerscales[4] = 1.0;
    optimizerscales[5] = 1.0;
    optimizerscales[6] = 1.0;
    optimizerscales[7] = 1.0;
    optimizerscales[8] = 1.0;
    optimizerscales[9]  = translscale;
    optimizerscales[10] = translscale;
    optimizerscales[11] = translscale;
    optimizer->SetScales( optimizerscales );
    optimizer->SetMaximumStepLength( 0.2  );
    optimizer->SetMinimumStepLength( 0.001 );
	optimizer->SetRelaxationFactor(0.5);
    optimizer->SetNumberOfIterations( steps );
	optimizer->MinimizeOn();

    // Registration
    registration->SetOptimizer(optimizer);
    registration->SetTransform(affinetransf);
	registration->SetFixedImageRegion(fixreg);
    registration->SetInitialTransformParameters(affinetransf->GetParameters());
    try {
        registration->StartRegistration();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(affinetransf);
    }
    affinetransf->SetParameters(registration->GetLastTransformParameters());
    std::cout << "\t\\-- Affine transform completed" << std::endl;
    std::cout << std::endl;

    return(affinetransf);
}

DeformationField BrainRegistration::DemonsRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int steps) {
	std::cout << "\tBrainRegistration::DemonsRegistration" << std::endl;
    /* Init */
    DemonsRegistrationFilterType::Pointer registration = DemonsRegistrationFilterType::New();
	ObserverType::Pointer observer = ObserverType::New();
	DeformationField deformation = DeformationFieldType::New();

	// Parameters
	registration->SetFixedImage( fiximg );
	registration->SetMovingImage( movimg );
	registration->SetNumberOfIterations( steps );
	registration->SetStandardDeviations( 1.0 );
	//registration->AddObserver( itk::IterationEvent(), observer );

	// Registration
	try {
		registration->Update();
	}
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(deformation);
    }

	deformation = registration->GetOutput();

	return(deformation);
}

RigidTransform BrainRegistration::MultiRigidRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels, int steps) {
	std::cout << "\tBrainRegistration::MultiRigidRegistration" << std::endl;
    /* Init */
    // Definitions
    const double translscale = 1.0 / 1000.0;
    ProbabilityRegionType fixreg = movimg->GetLargestPossibleRegion();
    // Observer
    ObserverType::Pointer observer = ObserverType::New();
    // Metric parameters
    PMattesMIMetricType::Pointer mattesmi_metric = PMattesMIMetricType::New();
	mattesmi_metric->SetNumberOfHistogramBins(50);
	mattesmi_metric->ComputeGradientOff();
	try {
        mattesmi_metric->SetNumberOfSpatialSamples(static_cast<int>(fixreg.GetNumberOfPixels()*0.5));
    }
    catch( itk::ExceptionObject ) {}
    // Transformations
	VersorType rotation;
	VectorType axis;
    RigidTransformType::Pointer rigidtransf = RigidTransformType::New();
    // Optimizers
    VersorOptimizerType::Pointer versor_optimizer = VersorOptimizerType::New();
    versor_optimizer->AddObserver( itk::IterationEvent(), observer );
    // Interpolator
    PInterpolatorType::Pointer interpolator = PInterpolatorType::New();
    // Registration framework
    PMultiResolutionRegistrationType::Pointer multires_registration = PMultiResolutionRegistrationType::New();
    multires_registration->SetMovingImage(movimg);
    multires_registration->SetFixedImage(fiximg);
    multires_registration->SetMetric(mattesmi_metric);
    multires_registration->SetInterpolator(interpolator);
    PMultiResolutionPyramidType::Pointer multires_pyramidfix = PMultiResolutionPyramidType::New();
    multires_registration->SetFixedImagePyramid(multires_pyramidfix);
    PMultiResolutionPyramidType::Pointer multires_pyramidmov = PMultiResolutionPyramidType::New();
    multires_registration->SetMovingImagePyramid(multires_pyramidmov);
    multires_registration->AddObserver( itk::IterationEvent(), observer );

    /* Rigid */
    std::cout << std::endl;
    std::cout << "\t\\-- Initializing rigid transform" << std::endl;
    // Registration initialization
	axis[0] = 0.0;
	axis[1] = 0.0;
	axis[2] = 1.0;
	rotation.Set( axis, 0.0 );
	rigidtransf->SetRotation( rotation );
    PTransformInitializerType::Pointer initializer = PTransformInitializerType::New();
    initializer->SetFixedImage(fiximg);
    initializer->SetMovingImage(movimg);
    initializer->SetTransform(rigidtransf);
    initializer->MomentsOn();
	//initializer->GeometryOn();
    initializer->InitializeTransform();

    // Optimizer parameters
    VersorScalesType versoroptimizerscales(rigidtransf->GetNumberOfParameters());
    versoroptimizerscales[0] = 1.0;
    versoroptimizerscales[1] = 1.0;
    versoroptimizerscales[2] = 1.0;
    versoroptimizerscales[3] = translscale;
    versoroptimizerscales[4] = translscale;
    versoroptimizerscales[5] = translscale;
    versor_optimizer->SetScales( versoroptimizerscales );
	versor_optimizer->SetMaximumStepLength( 0.2000  );
    versor_optimizer->SetMinimumStepLength( 0.0001 );
	versor_optimizer->SetRelaxationFactor(0.5);
    versor_optimizer->SetNumberOfIterations( steps );
	versor_optimizer->MinimizeOn();

    // Registration
    multires_registration->SetOptimizer(versor_optimizer);
    multires_registration->SetTransform(rigidtransf);
    multires_registration->SetFixedImageRegion(fixreg);
    multires_registration->SetNumberOfLevels(levels);
    multires_registration->SetInitialTransformParameters(rigidtransf->GetParameters());
    std::cout << "\t\\-- Starting rigid transform" << std::endl;
    try {
        multires_registration->StartRegistration();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(rigidtransf);
    }
    rigidtransf->SetParameters(multires_registration->GetLastTransformParameters());
    std::cout << "\t\\-- Rigid transform completed" << std::endl;
    std::cout << std::endl;

	return(rigidtransf);
}

AffineTransform BrainRegistration::MultiAffineRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels, int steps) {
    /* Rigid */
    RigidTransformType::Pointer rigidtransf = MultiRigidRegistration(fiximg, movimg, levels, steps);

	std::cout << "\tBrainRegistration::MultiAffineRegistration" << std::endl;
	/* Init */
    // Definitions
    const double translscale = 1.0 / 1000.0;
    IntensityRegionType fixreg = fiximg->GetLargestPossibleRegion();
    // Observer
    ObserverType::Pointer observer = ObserverType::New();
    // Metric parameters
    PMattesMIMetricType::Pointer mattesmi_metric = PMattesMIMetricType::New();
	mattesmi_metric->ComputeGradientOff();
	mattesmi_metric->SetNumberOfHistogramBins(50);
	try {
        mattesmi_metric->SetNumberOfSpatialSamples(static_cast<int>(fixreg.GetNumberOfPixels()*0.5));
    }
    catch( itk::ExceptionObject ) {}
    // Transformations
    AffineTransformType::Pointer affinetransf = AffineTransformType::New();
    // Optimizers
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );
    // Interpolator
    PInterpolatorType::Pointer interpolator = PInterpolatorType::New();
    // Registration framework
    PMultiResolutionRegistrationType::Pointer multires_registration = PMultiResolutionRegistrationType::New();
    multires_registration->SetMovingImage(movimg);
    multires_registration->SetFixedImage(fiximg);
    multires_registration->SetMetric(mattesmi_metric);
    multires_registration->SetInterpolator(interpolator);
    PMultiResolutionPyramidType::Pointer multires_pyramidfix = PMultiResolutionPyramidType::New();
    multires_registration->SetFixedImagePyramid(multires_pyramidfix);
    PMultiResolutionPyramidType::Pointer multires_pyramidmov = PMultiResolutionPyramidType::New();
    multires_registration->SetMovingImagePyramid(multires_pyramidmov);
    multires_registration->AddObserver( itk::IterationEvent(), observer );

    /* Affine */
    std::cout << std::endl;
    std::cout << "\t\\-- Starting affine transform" << std::endl;
	std::cout << std::endl;
    // Transform parameters
    affinetransf->SetCenter(rigidtransf->GetCenter());
    affinetransf->SetTranslation(rigidtransf->GetTranslation());
    affinetransf->SetMatrix(rigidtransf->GetMatrix());

    // Optimizer parameters
    OptimizerScalesType optimizerscales(affinetransf->GetNumberOfParameters());
    optimizerscales[0] = 1.0;
    optimizerscales[1] = 1.0;
    optimizerscales[2] = 1.0;
    optimizerscales[3] = 1.0;
    optimizerscales[4] = 1.0;
    optimizerscales[5] = 1.0;
    optimizerscales[6] = 1.0;
    optimizerscales[7] = 1.0;
    optimizerscales[8] = 1.0;
    optimizerscales[9]  = translscale;
    optimizerscales[10] = translscale;
    optimizerscales[11] = translscale;
    optimizer->SetScales( optimizerscales );
    optimizer->SetMaximumStepLength( 0.2 );
    optimizer->SetMinimumStepLength( 0.0001 );
    optimizer->SetNumberOfIterations( steps );

    // Registration
    multires_registration->SetOptimizer(optimizer);
    multires_registration->SetTransform(affinetransf);
	multires_registration->SetFixedImageRegion(fixreg);
    multires_registration->SetNumberOfLevels(levels);
    multires_registration->SetInitialTransformParameters(affinetransf->GetParameters());
    try {
        //registration->StartRegistration();
        multires_registration->StartRegistration();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(affinetransf);
    }
    affinetransf->SetParameters(multires_registration->GetLastTransformParameters());
    std::cout << "\t\\-- Affine transform completed" << std::endl;
    std::cout << std::endl;

    return(affinetransf);
}

DeformableTransform BrainRegistration::BSplineMultiRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, AffineTransform affinetransf, int levels, unsigned int ngridnodesdim, int steps) {
	std::cout << "\tBrainRegistration::BSplineMultiRegistration" << std::endl;
	/* Init */
	std::cout << "\t\\- Init" << std::endl;
    // Definitions
	unsigned int nbsplineparam;
	RegionType::SizeType gridsizeimg;
    RegionType::SizeType gridbordersize;
    RegionType::SizeType totalgridsize;
	OriginType origin;
    SpacingType spacing;
    ProbabilityRegionType fixreg;
    ProbabilityDirectionType griddir;
	ProbabilitySizeType fixsize;
	RegionType bsplinereg;

	// Grid estimation
    std::cout << "\t\t Grid estimation" << std::endl;
	origin = fiximg->GetOrigin();
    spacing = fiximg->GetSpacing();
    griddir = fiximg->GetDirection();
	fixreg = fiximg->GetLargestPossibleRegion();
	fixsize = fixreg.GetSize();
    gridsizeimg.Fill(ngridnodesdim);
    gridbordersize.Fill(3);    // Border for spline order = 3 ( 1 lower, 2 upper )
    totalgridsize = gridsizeimg + gridbordersize;
    bsplinereg.SetSize(totalgridsize);
    for(unsigned int r=0; r<3; r++) {
        spacing[r] *= static_cast<double>(fixsize[r] - 1)  /
                    static_cast<double>(gridsizeimg[r] - 1);
    }
    OriginType gridorig = origin - griddir * spacing;

    // Observer
    ObserverType::Pointer observer = ObserverType::New();

    // Metric parameters
    std::cout << "\t\t Metric parameters" << std::endl;
    PMattesMIMetricType::Pointer mattesmi_metric = PMattesMIMetricType::New();
	/*mattesmi_metric->ComputeGradientOff();
    mattesmi_metric->SetUseCachingOfBSplineWeights(false);  
	mattesmi_metric->SetNumberOfHistogramBins(50);
	mattesmi_metric->UseAllPixelsOn();
    mattesmi_metric->SetUseExplicitPDFDerivatives(false);*/
    try {
        mattesmi_metric->SetNumberOfSpatialSamples(static_cast<int>(fixreg.GetNumberOfPixels()*0.5));
    }
    catch( itk::ExceptionObject ) {}

    // Transformation
    std::cout << "\t\t Transformation" << std::endl;
    DeformableTransformType::Pointer bsplinetransf = DeformableTransformType::New();
	bsplinetransf->SetGridSpacing(spacing);
    bsplinetransf->SetGridOrigin(gridorig);
    bsplinetransf->SetGridRegion(bsplinereg);
    bsplinetransf->SetGridDirection(griddir);
    bsplinetransf->SetBulkTransform(affinetransf);
    nbsplineparam = bsplinetransf->GetNumberOfParameters();
    ParametersType initdeftransfparam(nbsplineparam);
    initdeftransfparam.Fill( 0.0 );
    bsplinetransf->SetParameters(initdeftransfparam);

    // Optimizer
    std::cout << "\t\t Optimizer" << std::endl;
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->AddObserver( itk::IterationEvent(), observer );
	OptimizerScalesType optimizerscales = OptimizerScalesType(nbsplineparam);
    optimizerscales.Fill( 1.0 );
    optimizer->SetScales( optimizerscales );
    optimizer->SetMaximumStepLength( 10.0 );
    optimizer->SetMinimumStepLength(  0.01 );
    optimizer->SetNumberOfIterations( steps );

    // Interpolator
    PInterpolatorType::Pointer interpolator = PInterpolatorType::New();

    // Registration framework
    std::cout << "\t\t Registration" << std::endl;
    PMultiResolutionRegistrationType::Pointer multires_registration = PMultiResolutionRegistrationType::New();
	PMultiResolutionPyramidType::Pointer multires_pyramidfix = PMultiResolutionPyramidType::New();
    multires_registration->SetFixedImagePyramid(multires_pyramidfix);
    PMultiResolutionPyramidType::Pointer multires_pyramidmov = PMultiResolutionPyramidType::New();
    multires_registration->SetMovingImagePyramid(multires_pyramidmov);
    multires_registration->SetMovingImage(movimg);
    multires_registration->SetFixedImage(fiximg);
    multires_registration->SetMetric(mattesmi_metric);
	multires_registration->SetOptimizer(optimizer);
    multires_registration->SetInterpolator(interpolator);
    multires_registration->SetNumberOfLevels(levels);
    multires_registration->SetTransform(bsplinetransf);
	multires_registration->AddObserver( itk::IterationEvent(), observer );

    /* BSpline */
    std::cout << std::endl;
    std::cout << "\t\\- Starting bspline transform" << std::endl;
    // Registration
    try {
      multires_registration->StartRegistration();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(bsplinetransf);
    }
    bsplinetransf->SetParameters(multires_registration->GetLastTransformParameters());
    std::cout << "\t\\- Bspline transform completed" << std::endl;
    std::cout << std::endl;

	return(bsplinetransf);
}

DeformableTransform BrainRegistration::BSplineMultiRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels, unsigned int ngridnodesdim, int steps) {
	/* Affine */
	AffineTransformType::Pointer affinetransf = MultiAffineRegistration(fiximg, movimg);

	DeformableTransform bsplinetransf = BSplineMultiRegistration(fiximg, movimg, affinetransf, levels, ngridnodesdim, steps);
	return(bsplinetransf);
}

DeformationField BrainRegistration::MultiDemonsRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels, int steps) {
	std::cout << "\tBrainRegistration::MultiDemonsRegistration" << std::endl;
    /* Init */
	unsigned int *nIterations = new unsigned int[levels];
	for(int i=0; i<levels; i++) nIterations[i] = steps;
	MultiDemonsRegistrationFilterType::Pointer registration = MultiDemonsRegistrationFilterType::New();
    DemonsRegistrationFilterType::Pointer demons = DemonsRegistrationFilterType::New();
	ObserverType::Pointer observer = ObserverType::New();
	DeformationField deformation = DeformationFieldType::New();

	// Parameters
	demons->SetStandardDeviations( 1.0 );
	registration->SetFixedImage( fiximg );
	registration->SetMovingImage( movimg );
	registration->SetNumberOfIterations( nIterations );
	registration->SetRegistrationFilter( demons );
	registration->SetNumberOfLevels( levels  );
	//registration->AddObserver( itk::IterationEvent(), observer );
	// Registration
	try {
		registration->Update();
	}
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(deformation);
    }

	deformation = registration->GetOutput();

	delete [] nIterations;

	return(deformation);
}

ProbabilityImage BrainRegistration::Subtraction(ProbabilityImage minuend, ProbabilityImage subtrahend, RigidTransform transformation) {
	std::cout << "\tBrainRegistration::Subtraction" << std::endl;
	// Init
	RigidTransformType::OutputVectorType half_translation;
	RigidTransformType::OutputVectorType inverse_half_translation;
	VersorType half_rotation;
	VersorType inverse_half_rotation;
	RigidTransformType::OutputVectorType::ValueType tx;
	RigidTransformType::OutputVectorType::ValueType ty;
	RigidTransformType::OutputVectorType::ValueType tz;
	VersorType::ValueType angle_radian;
	VersorType::VectorType versor_axis;
	double half_angle;
	double inverse_half_angle;


	// Transformation data
	std::cout << "\t\\- Getting the versor" << std::endl;
	VersorType versor =  transformation->GetVersor();
	AffineTransformType::OutputVectorType translation = transformation->GetTranslation();

	// Center
	std::cout << "\t\tCenter" << std::endl;
	BaseImageCalculatorType::Pointer m_BaseCalculator = BaseImageCalculatorType::New();
	m_BaseCalculator->SetImage(  minuend );		    
	m_BaseCalculator->Compute();
	BaseImageCalculatorType::VectorType baseCenter  = m_BaseCalculator->GetCenterOfGravity();

	AffineTransformType::CenterType centerFollowingImage;
	centerFollowingImage[0] = baseCenter[0];
	centerFollowingImage[1] = baseCenter[1];
	centerFollowingImage[2] = baseCenter[2];

	// Translation
	std::cout << "\t\tTranslation" << std::endl;
	tx = translation[0];
	ty = translation[1];
	tz = translation[2];

	half_translation[0] = tx/2;
	half_translation[1] = ty/2;
	half_translation[2] = tz/2;

	inverse_half_translation[0] = -half_translation[0];
	inverse_half_translation[1] = -half_translation[1];
	inverse_half_translation[2] = -half_translation[2];

	// Rotation
	std::cout << "\t\tRotation" << std::endl;
	angle_radian = versor.GetAngle();
	versor_axis = versor.GetAxis();
	VersorType::VectorType versor_axis_t;
	versor_axis_t[0] = versor_axis[0];
	versor_axis_t[1] = versor_axis[1];
	versor_axis_t[2] = versor_axis[2];
	std::cout << "\t\\--- Axis: " << versor_axis << std::endl;
	std::cout << "\t\\--- Angle: " << angle_radian << std::endl;
	std::cout << "\t\\--- Translation : " << half_translation << std::endl;
	half_angle = angle_radian/2;
	inverse_half_angle = -half_angle;
	

	half_rotation.Set(  versor_axis_t, half_angle  );
	inverse_half_rotation.Set( versor_axis_t,inverse_half_angle );		

	// Transforms
	std::cout << "\t\\- Transformations" << std::endl;
	RigidTransform halfRigid = RigidTransformType::New();
	halfRigid->SetCenter(transformation->GetCenter());
	halfRigid->SetRotation( half_rotation );
	halfRigid->SetTranslation( half_translation );
	
	RigidTransform inverseHalfRigid = RigidTransformType::New();
	inverseHalfRigid->SetCenter(transformation->GetCenter());
	inverseHalfRigid->SetRotation( inverse_half_rotation );
	inverseHalfRigid->SetTranslation(inverse_half_translation);

	// Halfway resample and subtraction
	std::cout << "\t\\- Halfway resample and subtraction" << std::endl;
	BSplineAtlasInterpolatorType::Pointer interpolator = BSplineAtlasInterpolatorType::New();

	ResampleAtlasFilterType::Pointer basalResample = ResampleAtlasFilterType::New();
	basalResample = ResampleAtlasFilterType::New();
	basalResample->SetInput(subtrahend);
	basalResample->SetSize(minuend->GetLargestPossibleRegion().GetSize());
	basalResample->SetOutputOrigin(minuend->GetOrigin());
	basalResample->SetOutputSpacing(minuend->GetSpacing());
	basalResample->SetOutputDirection(minuend->GetDirection());
	basalResample->SetInterpolator(interpolator);
	basalResample->SetDefaultPixelValue( 0 );
	basalResample->SetTransform(halfRigid);
	try {
		basalResample->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	ResampleAtlasFilterType::Pointer followResample = ResampleAtlasFilterType::New();
	followResample = ResampleAtlasFilterType::New();
	followResample->SetInput(minuend);
	followResample->SetSize(minuend->GetLargestPossibleRegion().GetSize());
	followResample->SetOutputOrigin(minuend->GetOrigin());
	followResample->SetOutputSpacing(minuend->GetSpacing());
	followResample->SetOutputDirection(minuend->GetDirection());
	followResample->SetInterpolator(interpolator);
	followResample->SetDefaultPixelValue( 0 );
	followResample->SetTransform(inverseHalfRigid);
	try {
		followResample->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	// We resample the subtraction
	std::cout << "\t\\- Subtraction resampling" << std::endl;
	SubtractionFilterType::Pointer subtractor = SubtractionFilterType::New();
	subtractor->SetInput1( followResample->GetOutput() );
	subtractor->SetInput2( basalResample->GetOutput() );
	subtractor->Update();

	ResampleAtlasFilterType::Pointer subResample = ResampleAtlasFilterType::New();
	subResample = ResampleAtlasFilterType::New();
	subResample->SetInput(subtractor->GetOutput());
	subResample->SetSize(minuend->GetLargestPossibleRegion().GetSize());
	subResample->SetOutputOrigin(minuend->GetOrigin());
	subResample->SetOutputSpacing(minuend->GetSpacing());
	subResample->SetOutputDirection(minuend->GetDirection());
	subResample->SetInterpolator(interpolator);
	subResample->SetDefaultPixelValue( 0 );
	subResample->SetTransform(halfRigid);
	try {
		subResample->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	return(subResample->GetOutput());
}

ProbabilityImage BrainRegistration::Warp(ProbabilityImage fiximg, ProbabilityImage movimg, DeformationField deformation) {

	WarperType::Pointer warper = WarperType::New(); 
	InterpolatorType::Pointer interpolator = InterpolatorType::New(); 
	
	warper->SetInput( fiximg ); 
	warper->SetInterpolator( interpolator ); 
	warper->SetOutputSpacing( movimg->GetSpacing() ); 
	warper->SetOutputOrigin( movimg->GetOrigin() );
	warper->SetOutputDirection( movimg->GetDirection() );
	warper->SetDeformationField( deformation );
	warper->Update();

	return( warper->GetOutput() );
}

MaskImage BrainRegistration::Warp(MaskImage fiximg, MaskImage movimg, DeformationField deformation) {

	WarperMaskType::Pointer warper = WarperMaskType::New(); 
	NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New(); 
	
	warper->SetInput( fiximg ); 
	warper->SetInterpolator( interpolator ); 
	warper->SetOutputSpacing( movimg->GetSpacing() ); 
	warper->SetOutputOrigin( movimg->GetOrigin() );
	warper->SetOutputDirection( movimg->GetDirection() );
	warper->SetDeformationField( deformation );
	warper->Update();

	return( warper->GetOutput() );
}

DeformationField BrainRegistration::Sum(std::vector<DeformationField> deformations) {
	std::cout << "\tBrainRegistration::Sum" << std::endl;
	VectorPixelType zero;
	zero.Fill(0);
	DeformationField average = DeformationFieldType::New();
    average->SetRegions( deformations[0]->GetLargestPossibleRegion() );
    average->SetSpacing( deformations[0]->GetSpacing() );
    average->SetOrigin( deformations[0]->GetOrigin() );
    average->SetDirection( deformations[0]->GetDirection() );
    average->Allocate();
    average->FillBuffer(zero);
    average->Update();
	std::vector<DeformationField>::iterator it;
	MultiplyConstantFilterType::Pointer constFilter = MultiplyConstantFilterType::New();
	
	it = deformations.begin();
	std::cout << "\t\\-- Sum of the images" << std::endl;
	for (it = deformations.begin(); it != deformations.end(); ++it) {
		AddFilterType:: Pointer add = AddFilterType::New();
		add->SetInput1(average);
        add->SetInput2(*it);
        add->Update();
		average = add->GetOutput();
	}
	
	return(average);
}

ProbabilityImage BrainRegistration::Jacobian(DeformationField deformation) {

	JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
	jacobianFilter->SetInput( deformation );
	jacobianFilter->Update();

	return( jacobianFilter->GetOutput() );
}

ProbabilityImage BrainRegistration::Divergence(DeformationField deformation) {

	/* Init */
	ProbabilityImage divergence = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	divergence->FillBuffer(0);
	ProbabilityIterator itDiv = ProbabilityIterator(divergence, divergence->GetLargestPossibleRegion());
	

	/* Vector decomposition */
	ProbabilityImage vx = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	vx->FillBuffer(0);
	ProbabilityImage vy = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	vy->FillBuffer(0);
	ProbabilityImage vz = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	vz->FillBuffer(0);

	DeformationIterator itDef = DeformationIterator(deformation, deformation->GetLargestPossibleRegion());
	ProbabilityIterator itVx = ProbabilityIterator(vx, vx->GetLargestPossibleRegion());
	ProbabilityIterator itVy = ProbabilityIterator(vy, vy->GetLargestPossibleRegion());
	ProbabilityIterator itVz = ProbabilityIterator(vz, vz->GetLargestPossibleRegion());

	for (itDef.GoToBegin(); !itDef.IsAtEnd(); ++itDef, ++itVx, ++itVy, ++itVz) {
		itVx.Set(itDef.Get()[0]);
		itVy.Set(itDef.Get()[1]);
		itVz.Set(itDef.Get()[2]);
	}

	/* Gradients */
	GradientFilterType::Pointer GradientVx = GradientFilterType::New();
	GradientVx->SetInput(vx);
	GradientVx->Update();
	GradientImage gx = GradientVx->GetOutput();
	GradientIterator itGx = GradientIterator(gx, gx->GetLargestPossibleRegion());
	GradientFilterType::Pointer GradientVy = GradientFilterType::New();
	GradientVy->SetInput(vy);
	GradientVy->Update();
	GradientImage gy = GradientVy->GetOutput();
	GradientIterator itGy = GradientIterator(gy, gy->GetLargestPossibleRegion());
	GradientFilterType::Pointer GradientVz = GradientFilterType::New();
	GradientVz->SetInput(vz);
	GradientVz->Update();
	GradientImage gz = GradientVz->GetOutput();
	GradientIterator itGz = GradientIterator(gz, gz->GetLargestPossibleRegion());

	for (itDiv.GoToBegin(); !itDiv.IsAtEnd(); ++itDiv, ++itGx, ++itGy, ++itGz)
        itDiv.Set(itGx.Get()[0]+itGy.Get()[1]+itGz.Get()[2]);

	return( divergence );
}
