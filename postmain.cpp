#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainregistration.h>

using namespace std;

int main(int argc, char **argv)
{
    std::cout << "/-------------------------------/" << std::endl;
    std::cout << "|            Posttool           |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;

    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baselineFolderName, baselineImagesFolderName, baselineTransformsFolderName, baselineSegmentationFolderName;
	FileNameType followupFolderName, followupImagesFolderName, followupTransformsFolderName, followupSegmentationFolderName;
	FileNameType subtractionFolderName, deformationFolderName;
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
            deformationFolderName = followupFolderName + "deformation/";
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
            deformationFolderName = followupFolderName + "deformation/";
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
            deformationFolderName = followupFolderName + "deformation/";
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
            deformationFolderName = followupFolderName + "deformation/";
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
            deformationFolderName = followupFolderName + "deformation/";
		case 8:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			baselineSegmentationFolderName = baselineFolderName + argv[5];
			followupSegmentationFolderName = followupFolderName + argv[5];
			subtractionFolderName = followupFolderName + argv[6];
			deformationFolderName = followupFolderName + argv[7];
		default:
            std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: PostTool baseline_folder followup_folder [images_folder] [transforms_folder] [segmentation_folder] [subtraction_folder] [deformation_folder]" << std::endl;
			return EXIT_FAILURE;
			break;
	}

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	//------> Deformation
	//-- Jacobian
	FileNameType jacobianName = deformationFolderName + "avg_multidemons_jacobian.nii.gz";
	ProbabilityImage jacobian = brainio->ReadProbabilityImage(jacobianName);
	FileNameType pdJacobianName = deformationFolderName + "pd_multidemons_jacobian.nii.gz";
	ProbabilityImage pdJacobian = brainio->ReadProbabilityImage(pdJacobianName);
	FileNameType t2JacobianName = deformationFolderName + "t2_multidemons_jacobian.nii.gz";
	ProbabilityImage t2Jacobian = brainio->ReadProbabilityImage(t2JacobianName);
	FileNameType flairJacobianName = deformationFolderName + "flair_multidemons_jacobian.nii.gz";
	ProbabilityImage flairJacobian = brainio->ReadProbabilityImage(flairJacobianName);
	//-- Divergence
	FileNameType divergenceName = deformationFolderName + "avg_multidemons_divergence.nii.gz";
	ProbabilityImage divergence = brainio->ReadProbabilityImage(divergenceName);
	FileNameType pdDivergenceName = deformationFolderName + "pd_multidemons_divergence.nii.gz";
	ProbabilityImage pdDivergence = brainio->ReadProbabilityImage(pdDivergenceName);
	FileNameType t2DivergenceName = deformationFolderName + "t2_multidemons_divergence.nii.gz";
	ProbabilityImage t2Divergence = brainio->ReadProbabilityImage(t2DivergenceName);
	FileNameType flairDivergenceName = deformationFolderName + "flair_multidemons_divergence.nii.gz";
	ProbabilityImage flairDivergence = brainio->ReadProbabilityImage(flairDivergenceName);

	//------> PD 
	//-- Baseline
	FileNameType pdBaselineName = followupImagesFolderName + "pd_moved.nii.gz";
	ProbabilityImage pdBaseline = brainio->ReadProbabilityImage(pdBaselineName);
	//-- FollowUp
	FileNameType pdFollowupName = followupImagesFolderName + "pd_corrected.nii.gz";
	ProbabilityImage pdFollowup = brainio->ReadProbabilityImage(pdFollowupName);
	//-- Positive activity
	FileNameType pdPositiveActivityName = subtractionFolderName + "pd_positive_activity.nii.gz";
	MaskImage pdPositiveActivity = brainio->ReadMaskImage(pdPositiveActivityName);

	//------> T2
	//-- Baseline
	FileNameType t2BaselineName = followupImagesFolderName + "t2_moved.nii.gz";
	ProbabilityImage t2Baseline = brainio->ReadProbabilityImage(t2BaselineName);
	//-- FollowUp
	FileNameType t2FollowupName = followupImagesFolderName + "t2_corrected.nii.gz";
	ProbabilityImage t2Followup = brainio->ReadProbabilityImage(t2FollowupName);
	//-- Positive activity
	FileNameType t2PositiveActivityName = subtractionFolderName + "t2_positive_activity.nii.gz";
	MaskImage t2PositiveActivity = brainio->ReadMaskImage(t2PositiveActivityName);

	//------> FLAIR
	//-- Baseline
	FileNameType flairBaselineName = followupImagesFolderName + "flair_moved.nii.gz";
	ProbabilityImage flairBaseline = brainio->ReadProbabilityImage(flairBaselineName);
	//-- FollowUp
	FileNameType flairFollowupName = followupImagesFolderName + "flair_registered.nii.gz";
	ProbabilityImage flairFollowup = brainio->ReadProbabilityImage(flairFollowupName);
	//-- Positive activity
	FileNameType flairPositiveActivityName = subtractionFolderName + "flair_positive_activity.nii.gz";
	MaskImage flairPositiveActivity = brainio->ReadMaskImage(flairPositiveActivityName);

    std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t          Joint mask           " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
    std::vector<MaskImage> masks;
    if (pdPositiveActivity != (MaskImage)NULL)
        masks.push_back(pdPositiveActivity);
    if (t2PositiveActivity != (MaskImage)NULL)
        masks.push_back(t2PositiveActivity);
    if (flairPositiveActivity != (MaskImage)NULL)
        masks.push_back(flairPositiveActivity);
    MaskImage intersection = BrainSegmentation::Intersection(masks);
    brainio->WriteMaskImage(
        subtractionFolderName + "joint_positive_activity.nii.gz",
        intersection
    );

    ConnectComponentFilterType::Pointer connectedFilter = ConnectComponentFilterType::New();
    connectedFilter->SetInput(intersection);
    connectedFilter->Update();

    ConnectedImage labels = connectedFilter->GetOutput();
	
    brainio->WriteConnectedImage(
        subtractionFolderName + "joint_labels_positive_activity.nii.gz",
        labels
    );

	
    std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t     Demons postprocessing     " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
    if (pdPositiveActivity != (MaskImage)NULL) {
        MaskImage pdLesions = BrainSegmentation::DeformationPostProcessing(pdPositiveActivity, pdJacobian, pdDivergence);
        brainio->WriteMaskImage(
            subtractionFolderName + "pd_demons_positive_activity.nii.gz",
            pdLesions
        );
    }
    if (t2PositiveActivity != (MaskImage)NULL) {
        MaskImage t2Lesions = BrainSegmentation::DeformationPostProcessing(t2PositiveActivity, t2Jacobian, t2Divergence);
        brainio->WriteMaskImage(
            subtractionFolderName + "t2_demons_positive_activity.nii.gz",
            t2Lesions
        );
    }
    if  (flairPositiveActivity != (MaskImage)NULL) {
        MaskImage flairLesions = BrainSegmentation::DeformationPostProcessing(flairPositiveActivity, flairJacobian, flairDivergence);
        brainio->WriteMaskImage(
            subtractionFolderName + "flair_demons_positive_activity.nii.gz",
            flairLesions
        );
    }
    MaskImage lesions = BrainSegmentation::DeformationPostProcessing(intersection, jacobian, divergence);
    brainio->WriteMaskImage(
        subtractionFolderName + "joint_demons_positive_activity.nii.gz",
        lesions
    );


	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t     Ratio postprocessing      " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
    masks.clear();
    if (pdPositiveActivity != (MaskImage)NULL) {
        std::cout << "\t/-- PD" << std::endl;
        MaskImage pdLesions = BrainSegmentation::OnurPostProcessing(pdPositiveActivity, pdBaseline, pdFollowup);
        std::cout << "\t\t/-- Saving" << std::endl;
		brainio->WriteMaskImage(
			subtractionFolderName + "pd_onur_positive_activity.nii.gz",
            pdLesions
		);
        std::cout << "\t\t/-- Pushing" << std::endl;
        masks.push_back(pdLesions);
	}
	if (t2PositiveActivity != (MaskImage)NULL) {
        std::cout << "\t/-- T2" << std::endl;
		MaskImage t2Lesions = BrainSegmentation::OnurPostProcessing(t2PositiveActivity, t2Baseline, t2Followup);
        std::cout << "\t\t/-- Saving" << std::endl;
		brainio->WriteMaskImage(
			subtractionFolderName + "t2_onur_positive_activity.nii.gz",
			t2Lesions
		);
        masks.push_back(t2Lesions);
	}
	if  (flairPositiveActivity != (MaskImage)NULL) {
        std::cout << "\t/-- FLAIR" << std::endl;
		MaskImage flairLesions = BrainSegmentation::OnurPostProcessing(flairPositiveActivity, flairBaseline, flairFollowup);
        std::cout << "\t\t/-- Saving" << std::endl;
		brainio->WriteMaskImage(
			subtractionFolderName + "flair_onur_positive_activity.nii.gz",
			flairLesions
		);
        masks.push_back(flairLesions);
	}
    std::cout << "\t/-- JOINT" << std::endl;
    MaskImage onurLesions = BrainSegmentation::Intersection(masks);
    std::cout << "\t\t/-- Saving" << std::endl;
	brainio->WriteMaskImage(
		subtractionFolderName + "joint_onur_positive_activity.nii.gz",
		onurLesions
	);

	delete(brainio);

    return EXIT_SUCCESS;

}
