#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainpreprocessing.h>
#include <brainregistration.h>
#include <brainsegmentation.h>

using namespace std;

int main(int argc, char **argv)
{
    std::cout << "/-------------------------------/" << std::endl;
    std::cout << "|            Subtool            |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;
    BrainIO *brainio = new BrainIO();


	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baselineFolderName, baselineImagesFolderName, baselineTransformsFolderName, baselineSegmentationFolderName;
	FileNameType followupFolderName, followupImagesFolderName, followupTransformsFolderName, followupSegmentationFolderName;
	FileNameType subtractionFolderName;
	switch (argc) {
		case 3:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
            baselineImagesFolderName = baselineFolderName + "preprocessed/";
            followupImagesFolderName = followupFolderName + "preprocessed/";
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            baselineSegmentationFolderName = baselineFolderName + "segmentation/";
            followupSegmentationFolderName = followupFolderName + "segmentation/";
            subtractionFolderName = followupFolderName + "subtraction/";
			break;
		case 4:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            baselineSegmentationFolderName = baselineFolderName + "segmentation/";
            followupSegmentationFolderName = followupFolderName + "segmentation/";
            subtractionFolderName = followupFolderName + "subtraction/";
			break;
		case 5:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
            baselineSegmentationFolderName = baselineFolderName + "segmentation/";
            followupSegmentationFolderName = followupFolderName + "segmentation/";
            subtractionFolderName = followupFolderName + "subtraction/";
		case 6:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			baselineSegmentationFolderName = baselineFolderName + argv[5];
			followupSegmentationFolderName = followupFolderName + argv[5];
            subtractionFolderName = followupFolderName + "subtraction/";
		case 7:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			baselineSegmentationFolderName = baselineFolderName + argv[5];
			followupSegmentationFolderName = followupFolderName + argv[5];
			subtractionFolderName = followupFolderName + argv[6];
		default:
			std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: SubtractionTool baseline_folder followup_folder [images_folder] [transforms_folder] [segmentation_folder] [subtraction_folder]" << std::endl; 
			return EXIT_FAILURE;
			break;
	}
    #ifdef WIN32
    CreateDirectory(subtractionFolderName.c_str(),NULL);
    #else
    mkdir(subtractionFolderName.c_str(), 0777);
    #endif


	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
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
	AffineTransform flairBaselineAffine;
	if (pdBaseline != (ProbabilityImage)NULL)
		flairBaselineAffine = brainio->ReadAffineTransform(flairBaselineTransformName);
	//-- FollowUp
	FileNameType flairFollowupName = followupImagesFolderName + "flair_corrected.nii.gz";
	ProbabilityImage flairFollowup = brainio->ReadProbabilityImage(flairFollowupName);
	FileNameType flairFollowupTransformName = followupTransformsFolderName + "affineFLAIRtoSpace.tfm";
	if (pdFollowup != (ProbabilityImage)NULL) {
		AffineTransform flairFollowupAffine = brainio->ReadAffineTransform(flairFollowupTransformName);
		flairFollowup = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairFollowup, flairFollowupAffine);
	}

	//------> Segmentation files
	//-- Baseline
	MaskImage baselineWMMask = brainio->ReadMaskImage(baselineSegmentationFolderName + "wm_mask.nii.gz");
	//-- FollowUp
	MaskImage followupWMMask = brainio->ReadMaskImage(followupSegmentationFolderName + "wm_mask.nii.gz");
	

	/*--- Registration + subtraction ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Registration          " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We set the fixed and moving images we just preprocessed
    std::cout << "\t/-- Affine registration" << std::endl;
	FileNameType affineName = baselineTransformsFolderName + "baseline2followup.tfm";
	AffineTransform affine = brainio->ReadAffineTransform(affineName);
	if (affine == (AffineTransform)NULL) {
		if ((pdFollowup != (ProbabilityImage)NULL) && (pdBaseline != (ProbabilityImage)NULL)) {
			affine = BrainRegistration::AffineRegistration(pdFollowup, pdBaseline, 3000);
			brainio->WriteAffineTransform(
				baselineTransformsFolderName + "baseline2followup.tfm",
				affine
			);
		} else if ((pdFollowup != (ProbabilityImage)NULL) && (flairBaseline != (ProbabilityImage)NULL)) {
			affine = BrainRegistration::AffineRegistration(pdFollowup, flairBaseline, 3000);
			brainio->WriteAffineTransform(
				baselineTransformsFolderName + "baseline2followup.tfm",
				affine
			);
		}
		else {
			std::cerr << "The aren't any images to perform subtraction" << std::endl; 
			return EXIT_FAILURE;
		}
	}

    std::cout << "\t/-- Image resampling" << std::endl;
	std::cout << "\t -- WM mask" << std::endl;
	MaskImage baselineWMMaskMoved = BrainRegistration::ResampleImage<AffineTransformType>(followupWMMask, baselineWMMask, affine);
	ProbabilityImage pdMoved;
	if ((pdFollowup != (ProbabilityImage)NULL) && (pdBaseline != (ProbabilityImage)NULL)) {
		pdMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, pdBaseline, affine);
		brainio->WriteProbabilityImage(
			followupImagesFolderName + "pd_moved.nii.gz",
			pdMoved
		);
	}
	std::cout << "\t -- T2" << std::endl;
	ProbabilityImage t2Moved;
	if (t2Baseline != (ProbabilityImage)NULL) {
		t2Moved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t2Baseline, affine);
		brainio->WriteProbabilityImage(
			followupImagesFolderName + "t2_moved.nii.gz",
			t2Moved
		);
	}
	std::cout << "\t -- FLAIR" << std::endl;
	ProbabilityImage flairMoved;
	if (flairBaseline != (ProbabilityImage)NULL) {
		if ((pdBaseline != (ProbabilityImage)NULL) && (pdFollowup != (ProbabilityImage)NULL)) {
            std::cout << "\t   /-- PD in baseline and followup" << std::endl;
			AffineTransform flairAffine = AffineTransformType::New();
			affine->Compose(flairBaselineAffine, true);
			flairMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairBaseline, affine);
			brainio->WriteProbabilityImage(
				followupImagesFolderName + "flair_moved.nii.gz",
				flairMoved
			);
		}
		else if (pdFollowup != (ProbabilityImage)NULL) {
            std::cout << "\t   /-- PD only in followup" << std::endl;
			flairMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairBaseline, affine);
			brainio->WriteProbabilityImage(
				followupImagesFolderName + "flair_moved.nii.gz",
				flairMoved
			);
		}
		else {
            std::cout << "\t   /-- No PD" << std::endl;
			flairMoved = BrainRegistration::ResampleImage<AffineTransformType>(flairFollowup, flairBaseline, affine);
			brainio->WriteProbabilityImage(
				followupImagesFolderName + "flair_moved.nii.gz",
				flairMoved
			);
		}
	}

	// Final ROI
	MaskImage wmMask = BrainSegmentation::Union(baselineWMMaskMoved, followupWMMask);
	brainio->WriteMaskImage(
		baselineSegmentationFolderName + "union_wm_mask.nii.gz",
		wmMask
	);

	// The last step is to apply the subtraction
    std::cout << "\t/-- Subtractions" << std::endl;
	ProbabilityImage subtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_subtraction.nii.gz");
	if ((subtractionPD == (ProbabilityImage)NULL) && (pdFollowup != (ProbabilityImage)NULL) && (pdMoved != (ProbabilityImage)NULL)) {
		subtractionPD = BrainSegmentation::Subtraction(pdFollowup, pdMoved);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_subtraction.nii.gz",
			subtractionPD
		);
	}

	ProbabilityImage subtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_subtraction.nii.gz");
	if ((subtractionT2 == (ProbabilityImage)NULL) && (t2Followup != (ProbabilityImage)NULL) && (t2Moved != (ProbabilityImage)NULL)) {
		subtractionT2 = BrainSegmentation::Subtraction(t2Followup, t2Moved);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_subtraction.nii.gz",
			subtractionT2
		);
	}

	ProbabilityImage subtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_subtraction.nii.gz");
	if ((subtractionFLAIR == (ProbabilityImage)NULL) && (flairFollowup != (ProbabilityImage)NULL) && (flairMoved != (ProbabilityImage)NULL)) {
		subtractionFLAIR = BrainSegmentation::Subtraction(flairFollowup, flairMoved);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_subtraction.nii.gz",
			subtractionFLAIR
		);
	}
	

	// Masking with brainmask
	ProbabilityImage maskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_masked_subtraction.nii.gz");
	if ((maskedSubtractionPD == (ProbabilityImage)NULL) && (subtractionPD != (ProbabilityImage)NULL)) {
		maskedSubtractionPD = BrainPreprocessing<ProbabilityImageType>::MaskImage(subtractionPD, wmMask);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_masked_subtraction.nii.gz",
			maskedSubtractionPD
		);
	}

	ProbabilityImage maskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_masked_subtraction.nii.gz");
	if ((maskedSubtractionT2 == (ProbabilityImage)NULL) && (subtractionT2 != (ProbabilityImage)NULL)) {
		maskedSubtractionT2 = BrainPreprocessing<ProbabilityImageType>::MaskImage(subtractionT2, wmMask);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_masked_subtraction.nii.gz",
			maskedSubtractionT2
		);
	}

	ProbabilityImage maskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_masked_subtraction.nii.gz");
	if ((maskedSubtractionFLAIR == (ProbabilityImage)NULL) && (subtractionFLAIR != (ProbabilityImage)NULL)) {
		maskedSubtractionFLAIR = BrainPreprocessing<ProbabilityImageType>::MaskImage(subtractionFLAIR, wmMask);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_masked_subtraction.nii.gz",
			maskedSubtractionFLAIR
		);
	}

	/*--- Lesion detection ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t       Lesion detection        " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// Subtraction smoothing
	ProbabilityImage gaussSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_smoothed_subtraction.nii.gz");
	if ((gaussSubtractionPD == (ProbabilityImage)NULL) && (maskedSubtractionPD != (ProbabilityImage)NULL)) {
		gaussSubtractionPD = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedSubtractionPD, 0.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_smoothed_subtraction.nii.gz",
			gaussSubtractionPD
		);
	}

	ProbabilityImage gaussSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_smoothed_subtraction.nii.gz");
	if ((gaussSubtractionT2 == (ProbabilityImage)NULL) && (maskedSubtractionT2 != (ProbabilityImage)NULL)) {
		ProbabilityImage gaussSubtractionT2 = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedSubtractionT2, 0.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_smoothed_subtraction.nii.gz",
			gaussSubtractionT2
		);
	}

	ProbabilityImage gaussSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_smoothed_subtraction.nii.gz");
	if ((gaussSubtractionFLAIR == (ProbabilityImage)NULL) && (maskedSubtractionFLAIR != (ProbabilityImage)NULL)) {
		ProbabilityImage gaussSubtractionFLAIR = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedSubtractionFLAIR, 0.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_smoothed_subtraction.nii.gz",
			gaussSubtractionFLAIR
		);
	}


	// Simple thresholding of the smoothed subtractions
	MaskImage positiveActivityPD = brainio->ReadMaskImage(subtractionFolderName + "pd_positive_activity.nii.gz");
	if ((positiveActivityPD == (MaskImage)NULL) && (gaussSubtractionPD != (ProbabilityImage)NULL)) {
		MaskImage positiveActivityPD = BrainSegmentation::GetPositiveActivity(gaussSubtractionPD);
		brainio->WriteMaskImage(
			subtractionFolderName + "pd_positive_activity.nii.gz",
			positiveActivityPD
		);
	}

	MaskImage positiveActivityT2 = brainio->ReadMaskImage(subtractionFolderName + "t2_positive_activity.nii.gz");
	if ((positiveActivityT2 == (MaskImage)NULL) && (gaussSubtractionT2 != (ProbabilityImage)NULL)) {
		MaskImage positiveActivityT2 = BrainSegmentation::GetPositiveActivity(gaussSubtractionT2);
		brainio->WriteMaskImage(
			subtractionFolderName + "t2_positive_activity.nii.gz",
			positiveActivityT2
		);
	}

	MaskImage positiveActivityFLAIR = brainio->ReadMaskImage(subtractionFolderName + "flair_positive_activity.nii.gz");
	if ((positiveActivityFLAIR == (MaskImage)NULL) && (gaussSubtractionFLAIR != (ProbabilityImage)NULL)) {
		MaskImage positiveActivityFLAIR = BrainSegmentation::GetPositiveActivity(gaussSubtractionFLAIR);
		brainio->WriteMaskImage(
			subtractionFolderName + "flair_positive_activity.nii.gz",
			positiveActivityFLAIR
		);
	}


	delete(brainio);

    return EXIT_SUCCESS;

}
