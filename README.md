# BrainTools: C++ toolbox based on ITK for brain MRI processing.

__BrainTools__ is a C++ repository with several fully automated tools to process brain MRI data from MS patients. Currently it implements _bias correction_ using the __N4__ algorithm; _atlas registration_ using affine and b-splines transformations; _tissue segmentation_ based on <a name="deref1">[[1](#ref1)]</a>; _lesion segmentation and detection_ also based on [[1](#ref1)]; and, _longitudinal lesion analysis_ based on <a name="deref2">[[2](#ref2)]</a><a name="deref3">[[3](#ref3)]</a>.

# Dependencies:
+ This code supports version **3.20** of ITK. It has not been coded yet for newer versions and it won't compile with them. Also, the CMakelists might need some tweaking to define the ITK path depending on your OS.
+ CMake is needed to create a Makefile or project for your desired IDE. We recommend using CMake version 3+ for better compatibility.
+ In order to complement the preprocessing steps with skull stripping we recommend installing a separate toolbox. For instance, [BET](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET) is a widely used toolbox in medical imaging. In the example Python script we provide in this README we decided to use [ROBEX](https://sites.google.com/site/jeiglesias/ROBEX) due to their improved results and better integration with our pipeline.
+ If you decide to use the example script, you will also need python 2.7.
+ The **ICBM 452** atlas is not provided with this software. You can dowload it from the [LONI website](http://www.loni.usc.edu/atlases/Atlas_Detail.php?atlas_id=6). Once downloaded, the template file and the probabilistic maps should be located in the directory where you have your compiled files.

# Defaults
To use these tools you will need to acquire T1-w, T2-w, PD-w and T2-FLAIR-w images for each timepoint and patient. In order to automatically detect each image we assume that the image name will contain any of the substrings defined in the following table:

Brain mask |  T1-w  | T2-w | PD-w | T2-FLAIR-w
---------- | ------ | ---- | ---- | ----------
brain_mask | MPRAGE |  T2  |  DP  | flair
brainmask  | mprage |  t2  |  PD  | FLAIR
BrainMask  |  MPR   |      |  dp  | dark_fluid
           |  mpr   |      |  pd  | darkfluid
           |   T1   |      |      |
           |   t1   |      |      |

We have also defined default names for the following variables used when executing the tools:
+ **images_folder**: preprocessed/
+ **transforms_folder**: transforms/
+ **atlas_folder**: atlas/
+ **segmentation_folder**: segmentation/
+ **subtraction_folder**: subtraction/
+ **deformation_folder**: deformation/

# Bias correction
The name of the tool that implements this step is PreTool. The syntax for this tool is as follows:

```bash
PreTool baseline_folder followup_folder
```

This preprocessing tool is adapted to longitudinal analysis, thus, both followup and baseline folder must be specified. Also, we assume that a brainmask is located inside the patient's folder.
The output will be written to the _preprocessed/_ subdirectory inside the timepoint's folder.

To adapt that tool check the file check the file __premain.cpp__ and make the necessary adjustments.

# Coregistration
The name of the tool that implements this step is CoregTool. The syntax for this tool is as follows:

```bash
CoregTool base_folder [images_folder] [transformations_folder]
```

This tool coregisters the preprocessed FLAIR and T1 images inside *base_folder/images_folder*. The affine transformations are then stored inside *base_folder/transforms_folder*.

To adapt that tool check the file __coregmain.cpp__ and make the necessary adjustments.

# Atlas registration
The name of the tool that implements this step is AtlasTool. The syntax for this tool is as follows:

```bash
AtlasTool base_folder [images_folder] [transforms_folder] [atlas_folder]
```

This tool registers the **ICBM 452** atlas to a preprocessed T1 image inside *base_folder/images_folder* with its coregistration transformations inside *base_folder/transforms_folder*. As such, you will need to run PreTool and CoregTool first. Also, we assume that a brainmask is located inside the patient's folder.

The registered probabilistic maps will be written to the *base_folder/atlas_folder* subdirectory inside the timepoint's folder.

To adapt that tool check the file __atlasmain.cpp__ and make the necessary adjustments.

# Tissue and lesion segmentation
The name of the tool that implements this step is TissueTool. The syntax for this tool is as follows:

```bash
TissueTool base_folder [images_folder] [transforms_folder] [atlas_folder] [segmentation_folder]
```

This tool segments tissues using T1-w, T2-w, PD-w images stored in *base_folder/images_folder* with its coregistration transformations inside *base_folder/transforms_folder*. The segmentation based on [[1](#ref1)] uses an atlas located inside *base_folder/atlas_folder* to create an initial tissue segmentation with lesion voxels missclassified as tissue. This segmentation also includes a *partial volume* class that should be reclassified as either *grey matter* or *cerebro-spinal fluid*. After this initial segmentation, lesions are segmented using the T2-FLAIR-w image.

Finally, the tissue and lesion masks are stored inside *base_folder/segmentation_folder*.

To adapt that tool check the file check the file __tissuemain.cpp__ and make the necessary adjustments.

# Subtraction
The name of the tool that implements this step is SubtractionTool. The syntax for this tool is as follows:

```bash
SubtractionTool baseline_folder followup_folder [images_folder] [transforms_folder] [segmentation_folder] [subtraction_folder]
```

This tool applies a subtraction operation between the PD-w, T2-w and T2-FLAIR-w images inside the *baseline_folder/images_folder* and *followup_folder/images_folder* after applying and affine registration. Afterwards, a gaussian smoothing is applied and an initial new lesions mask is computed applying a threshold. All results are then saved inside *followup_folder/subtraction_folder*.

To adapt that tool check the file check the file __submain.cpp__ and make the necessary adjustments.

# Demons registration and deformation analysis

The name of the tool that implements this step is DeformableTool. The syntax for this tool is as follows:

```bash
DeformableTool baseline_folder followup_folder [images_folder] [transforms_folder] [deformation_folder]
```

This tool applies a multiresolution Demons algorithm from ITK between the PD-w, T2-w and T2-FLAIR-w images inside the *baseline_folder/images_folder* and *followup_folder/images_folder* after obtaining and initial affine transformation between them. For each deformation field, divergence and jacobian images are computed. All these outputs are stored in *followup_folder/deformation_folder*.

To adapt that tool check the file check the file __defomain.cpp__ and make the necessary adjustments.

# Final new lesion segmentation

The name of the tool that implements this step is PostTool. The syntax for this tool is as follows:

```bash
PostTool baseline_folder followup_folder [images_folder] [transforms_folder] [segmentation_folder] [subtraction_folder] [deformation_folder]
```

This tool obtains the final mask for the new and enhancing lesions using the initial mask stored in *followup_folder/subtraction_folder*, the deformation analysis from *followup_folder/deformation_folder* and the preprocessed images from *baseline_folder/images_folder* and *followup_folder/images_folder*. One of the final masks is computed using [[2](#ref2)], while the other is obtained using [[3](#ref3)]
The resulting masks are saved in *followup_folder/subtraction_folder*.

To adapt that tool check the file check the file _postmain.cpp__ and make the necessary adjustments.

# Example Python script
This script fully integrates the pipeline described in [[2](#ref2)] using ROBEX with an Ubuntu 14.04 machine. We assume that each patient has a folder with a numerical code inside /DATA/ and for each patient's folder each timepoint is labeled time__n__ where __n__ is the timepoint. This way, we can analyse differences between different timepoints. Terminal output is redirected to  _log.txt_ and any error is redirected to the file _error.log_ inside the patient's folder.


```python
import os
import sys
import getopt
from subprocess import call
import re


def main():
    folder_name = '/DATA/'

    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ('-h', '--help'):
            print('Help text')
            sys.exit(0)
        elif o in ('-f', '--folder'):
            folder_name = a

    patients = sorted(os.listdir(folder_name))
    for patient in patients:
        patient_folder = '%s%s/' % (folder_name, patient)
        print ''
        print '<---- Processing %s ---->' % (patient)
        print ''

        log_file = '%slog.txt' % (patient_folder)
        error_file = '%serror.txt' % (patient_folder)
        timepoints = [timepoint for timepoint in sorted(os.listdir(patient_folder)) if re.search(r'time(\d)*',timepoint)]
        with open(log_file, 'w') as f, open(error_file, 'w') as err:
           for timepoint in timepoints:
              is_timepoint = re.search(r'(time)(\d)*',timepoint)
              if is_timepoint:
                  timepoint_folder = '%s%s/' % (patient_folder, timepoint)
                  files = sorted(os.listdir(timepoint_folder))
                  mask_in_folder = False 
                    for file in files:
                    mask = re.search(r'(\w|\W)*(brainmask)(\w|\W)*',file)
                        if mask:
                            mask_in_folder = True
                        if not mask_in_folder:
                            for file in files:
                                pd_in_folder = re.search(r'(\w|\W)*(pd|dp|PD|DP)(\w|\W)*',file)
                                if pd_in_folder:
                                    pd_image = '%s%s' % (timepoint_folder, file)
                                    pd_image_stripped = '%sstripped.nii' % (timepoint_folder)
                                    mask_image = '%sbrainmask.nii' % (timepoint_folder)
                                    command = './ROBEX %s %s %s' % (pd_image, pd_image_stripped, mask_image)
                                    call([
                                        command
                                    ], stdout=f, stderr=err, shell=True)
                                    os.remove(pd_image_stripped)
             for baseline, followup in zip(timepoints[:-1], timepoints[1:]):
                 baseline_folder = '%s%s/' % (patient_folder, baseline)
                 followup_folder = '%s%s/' % (patient_folder, followup)
                 command = './PreTool %s %s' % (baseline_folder, followup_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './CoregTool %s' % (baseline_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './CoregTool %s' % (followup_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './AtlasTool %s' % (baseline_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './AtlasTool %s' % (followup_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './TissueTool %s' % (baseline_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './TissueTool %s' % (followup_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './SubtractionTool %s %s' % (baseline_folder, followup_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './DeformableTool %s %s' % (baseline_folder, followup_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
                 command = './PostTool %s %s' % (baseline_folder, followup_folder)
                 call([
                     command
                 ], stdout=f, stderr=err, shell=True)
if __name__ == '__main__':
    main()
```

<a name="ref1">[[1](#deref1)]</a> M. Cabezas, A. Oliver, E. Roura, J. Freixenet, J.C Vilanova, Ll. Ramió-Torrentà, A. Rovira, X. Lladó. [_Automatic multiple sclerosis lesion detection in brain MRI by FLAIR thresholding_](http://dx.doi.org/). __Computer Methods and Programs in Biomedicine__, 115(3), pp. 147-161. 2014

<a name="ref2">[[2](#deref2)]</a> M. Cabezas, J.F. Corral, A. Oliver, Y. Diez, M. Tintore, C. Auger, X. Montalban, X. Lladó, D. Pareto, A. Rovira. [_Automatic multiple sclerosis lesion detection in brain MRI by FLAIR thresholding_](http://dx.doi.org/10.1016/j.cmpb.2014.04.006). __American Journal of Neuroradiology__, to appear. 2016 

<a name="ref3">[[3](#deref3)]</a> O. Ganiler, A. Oliver, Y. Díez, J. Freixenet, J.C. Vilanova, B.Beltrán, Ll. Ramió-Torrentà, A. Rovira, and X. Lladó. [_A subtraction pipeline for automatic detection of new appearing multiple sclerosis lesions in longitudinal studies_](http://dx.doi.org/10.1007/s00234-014-1343-1). __Neuroradiology__, 56(5), pp. 363-374. 2014
