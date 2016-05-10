#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainpreprocessing.h>
#include <brainsegmentation.h>
#include <brainregistration.h>

using namespace std;

int main(int argc, char **argv)
{
    std::cout << "/-------------------------------/" << std::endl;
    std::cout << "|            Defotool           |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;

    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baselineFolderName, baselineImagesFolderName, baselineTransformsFolderName;
	FileNameType followupFolderName, followupImagesFolderName, followupTransformsFolderName;
	FileNameType deformationFolderName;

	switch (argc) {
		case 3:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
            baselineImagesFolderName = baselineFolderName + "preprocessed/";
            followupImagesFolderName = followupFolderName + "preprocessed/";
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            deformationFolderName = followupFolderName + "deformation/";
			break;
		case 4:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            deformationFolderName = followupFolderName + "deformation/";
			break;
		case 5:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
            deformationFolderName = followupFolderName + "deformation/";
		case 6:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			deformationFolderName = followupFolderName + argv[5];
		default:
			std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: SubtractionTool baseline_folder followup_folder [images_folder] [transforms_folder] [deformation_folder]" << std::endl; 
			return EXIT_FAILURE;
			break;
	}
    #ifdef WIN32
    CreateDirectory(deformationFolderName.c_str(),NULL);
    #else
    mkdir(deformationFolderName.c_str(), 0777);
    #endif
	FileNameType baselineSearchName = baselineFolderName + "*.nii";
	FileNameType followupSearchName = followupFolderName + "*.nii";
	std::vector<FileNameType> toFindVector;

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	//------> BRAINMASK
	toFindVector.clear();
	toFindVector.push_back("brain_mask");
	toFindVector.push_back("brainmask");
	toFindVector.push_back("BrainMask");
	//-- Baseline
	FileNameType maskBaselineName = BrainIO::SearchForFile(baselineSearchName, toFindVector);
	if (!maskBaselineName.empty())
		maskBaselineName = baselineFolderName + maskBaselineName;
	else
	{
		std::cerr << "ERROR: There is no brain mask with any of these names in folder "
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
		return EXIT_FAILURE;
	}
	MaskImage maskBaseline = brainio->ReadMaskImage(maskBaselineName);
	//-- FollowUp
	FileNameType maskFollowupName = BrainIO::SearchForFile(followupSearchName, toFindVector);
	if (!maskFollowupName.empty())
		maskFollowupName = followupFolderName + maskFollowupName;
	else
	{
		std::cerr << "ERROR: There is no brain mask with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
		return -1;
	}
	MaskImage maskFollowup = brainio->ReadMaskImage(maskFollowupName);

	//------> PD 
	//-- Baseline
	FileNameType pdBaselineName = baselineImagesFolderName + "pd_corrected_matched.nii.gz";
	ProbabilityImage pdBaseline = brainio->ReadProbabilityImage(pdBaselineName);
	//-- FollowUp
	FileNameType pdFollowupName = followupImagesFolderName + "pd_corrected.nii.gz";
	ProbabilityImage pdFollowup = brainio->ReadProbabilityImage(pdFollowupName);

	//------> T2
	//-- Baseline
	FileNameType t2BaselineName = baselineImagesFolderName + "t2_corrected_matched.nii.gz";
	ProbabilityImage t2Baseline = brainio->ReadProbabilityImage(t2BaselineName);
	//-- FollowUp
	FileNameType t2FollowupName = followupImagesFolderName + "t2_corrected.nii.gz";
	ProbabilityImage t2Followup = brainio->ReadProbabilityImage(t2FollowupName);

	//------> FLAIR
	//-- Baseline
	FileNameType flairBaselineName = baselineImagesFolderName + "flair_corrected_matched.nii.gz";
	ProbabilityImage flairBaseline = brainio->ReadProbabilityImage(flairBaselineName);
	FileNameType flairBaselineTransformName = baselineTransformsFolderName + "affineFLAIRtoSpace.tfm";
	AffineTransform flairBaselineAffine = brainio->ReadAffineTransform(flairBaselineTransformName);
	//-- FollowUp
	FileNameType flairFollowupName = followupImagesFolderName + "flair_corrected.nii.gz";
	ProbabilityImage flairFollowup = brainio->ReadProbabilityImage(flairFollowupName);
	FileNameType flairFollowupTransformName = followupTransformsFolderName + "affineFLAIRtoSpace.tfm";
	AffineTransform flairFollowupAffine = brainio->ReadAffineTransform(flairFollowupTransformName);
	if (flairFollowupAffine != (AffineTransform)NULL)
		flairFollowup = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairFollowup, flairFollowupAffine);

	//------> Affine from baseline to follow-up
	FileNameType affineName = baselineTransformsFolderName + "baseline2followup.tfm";
	AffineTransform affine2followup = brainio->ReadAffineTransform(affineName);
	
	/*--- Moving the baseline images ---*/
	MaskImage maskMoved = BrainRegistration::ResampleImage<AffineTransformType>(maskFollowup, maskBaseline, affine2followup);
	ProbabilityImage pdMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, pdBaseline, affine2followup);
	ProbabilityImage t2Moved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t2Baseline, affine2followup);
	if (flairBaselineAffine != (AffineTransform)NULL) {
		AffineTransformType::MatrixType flairBaselineTransformMatrix = flairBaselineAffine->GetMatrix();
		AffineTransformType::MatrixType baseline2FolloupMatrix = affine2followup->GetMatrix();
		AffineTransformType::MatrixType flairAffineTransformMatrix = baseline2FolloupMatrix * flairBaselineTransformMatrix;
		AffineTransform flairAffine = AffineTransformType::New();
		affine2followup->SetMatrix(flairAffineTransformMatrix);
	}
	ProbabilityImage flairMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairBaseline, affine2followup);
	
	/*--- Masking all images ---*/
	MaskImage mask = BrainSegmentation::Intersection(maskFollowup, maskMoved);
	ProbabilityImage pdMaskedBaseline = BrainPreprocessing<ProbabilityImageType>::MaskImage(pdMoved, mask);
	ProbabilityImage t2MaskedBaseline = BrainPreprocessing<ProbabilityImageType>::MaskImage(t2Moved, mask);
	ProbabilityImage flairMaskedBaseline = BrainPreprocessing<ProbabilityImageType>::MaskImage(flairMoved, mask);

	ProbabilityImage pdMaskedFollowup = BrainPreprocessing<ProbabilityImageType>::MaskImage(pdFollowup, mask);
	ProbabilityImage t2MaskedFollowup = BrainPreprocessing<ProbabilityImageType>::MaskImage(t2Followup, mask);
	ProbabilityImage flairMaskedFollowup = BrainPreprocessing<ProbabilityImageType>::MaskImage(flairFollowup, mask);

	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t   Demons multi Registration   " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;

	std::cout << "\tRegistration" << std::endl;
	DeformationField pdMRDeformation = brainio->ReadDeformationField(deformationFolderName + "pd_multidemons_deformation.nii.gz");
	if (pdMRDeformation == (DeformationField)NULL) {
        std::cout << "\t/-- Processing PD" << std::endl;
		pdMRDeformation = BrainRegistration::MultiDemonsRegistration(pdMaskedFollowup, pdMaskedBaseline);
        std::cout << "\t/---- Writing the results" << std::endl;
		brainio->WriteDeformationField(
			deformationFolderName + "pd_multidemons_deformation.nii.gz",
			pdMRDeformation
		);
	}

	DeformationField t2MRDeformation = brainio->ReadDeformationField(deformationFolderName + "t2_multidemons_deformation.nii.gz");
	if (t2MRDeformation == (DeformationField)NULL) {
        std::cout << "\t/-- Processing T2" << std::endl;
		t2MRDeformation = BrainRegistration::MultiDemonsRegistration(t2MaskedFollowup, t2MaskedBaseline);
        std::cout << "\t/---- Writing the results" << std::endl;
		brainio->WriteDeformationField(
			deformationFolderName + "t2_multidemons_deformation.nii.gz",
			t2MRDeformation
		);
	}

	DeformationField flairMRDeformation = brainio->ReadDeformationField(deformationFolderName + "flair_multidemons_deformation.nii.gz");
	if (flairMRDeformation == (DeformationField)NULL) {
        std::cout << "\t/-- Processing FLAIR" << std::endl;
		flairMRDeformation = BrainRegistration::MultiDemonsRegistration(flairMaskedFollowup, flairMaskedBaseline);
        std::cout << "\t/---- Writing the results" << std::endl;
		brainio->WriteDeformationField(
			deformationFolderName + "flair_multidemons_deformation.nii.gz",
			flairMRDeformation
		);
	}

	std::vector<DeformationField> deformations;
	if (pdMRDeformation != (DeformationField)NULL)
		deformations.push_back(pdMRDeformation);
	if (t2MRDeformation != (DeformationField)NULL)
		deformations.push_back(t2MRDeformation);
	if (flairMRDeformation != (DeformationField)NULL)
		deformations.push_back(flairMRDeformation);
	DeformationField avgMRDeformation = BrainRegistration::Sum(deformations);
	brainio->WriteDeformationField(
		deformationFolderName + "avg_multidemons_deformation.nii.gz",
		avgMRDeformation
	);

	std::cout << "\tJacobian" << std::endl;
	if (pdMRDeformation != (DeformationField)NULL) {
		ProbabilityImage pdJacobian = BrainRegistration::Jacobian(pdMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "pd_multidemons_jacobian.nii.gz",
			pdJacobian
		);
	}

	if (t2MRDeformation != (DeformationField)NULL) {
		ProbabilityImage t2Jacobian = BrainRegistration::Jacobian(t2MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t2_multidemons_jacobian.nii.gz",
			t2Jacobian
		);
	}

	if (flairMRDeformation != (DeformationField)NULL) {
		ProbabilityImage flairJacobian = BrainRegistration::Jacobian(flairMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "flair_multidemons_jacobian.nii.gz",
			flairJacobian
		);
	}

	ProbabilityImage avgJacobian = BrainRegistration::Jacobian(avgMRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg_multidemons_jacobian.nii.gz",
		avgJacobian
	);

	std::cout << "\tDivergence" << std::endl;
	if (pdMRDeformation != (DeformationField)NULL) {
		ProbabilityImage pdDivergence = BrainRegistration::Divergence(pdMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "pd_multidemons_divergence.nii.gz",
			pdDivergence
		);
	}

	if (t2MRDeformation != (DeformationField)NULL) {
		ProbabilityImage t2Divergence = BrainRegistration::Divergence(t2MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t2_multidemons_divergence.nii.gz",
			t2Divergence
		);
	}

	if (flairMRDeformation != (DeformationField)NULL) {
		ProbabilityImage flairDivergence = BrainRegistration::Divergence(flairMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "flair_multidemons_divergence.nii.gz",
			flairDivergence
		);
	}

	ProbabilityImage avgDivergence = BrainRegistration::Divergence(avgMRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg_multidemons_divergence.nii.gz",
		avgDivergence
	);

	delete(brainio);

    return EXIT_SUCCESS;

}
