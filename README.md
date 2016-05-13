# BrainTools: C++ code based on ITK for brain MRI processing.

__BrainTools__ is a C++ repository with several fully automated tools to process brain MRI data from MS patients. Currently it implements _bias correction_ using the __N4__ algorithm; _atlas registration_ using affine and b-splines transformations; _tissue segmentation_ based on [[1](#ref1)]; _lesion segmentation and detection_ also based on [[1](#ref1)]; and, _longitudinal lesion analysis_ based on [[2](#ref2)].

# Dependencies:
+ This code supports version **3.20** of ITK. It has not been coded yet for newer versions and it won't compile with them. Also, the CMakelists might need some tweaking to define the ITK path depending on your OS.
+ CMake is needed to create a Makefile or project for your desired IDE. We recommend using CMake version 3+ for better compatibility.
+ In order to complement the preprocessing steps with skull stripping we recommend installing a separate toolbox. For instance, [BET](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET) is a widely used toolbox in medical imaging. In the example Python script we provide in this README we decided to use [ROBEX](https://sites.google.com/site/jeiglesias/ROBEX) due to their improved results and better integration with our pipeline.
+ If you decide to use the example script, you will also need python 2.7.
 


# Example Python script
This script fully integrates the pipeline described in [[2](#ref2)] using ROBEX with an Ubuntu 14.04 machine. We assume that each patient has folder with a numerical code and inside this folder each timepoint is labeled time__n__ where __n__ is the timepoint. This way, we can analyse differences between different timepoints. Terminal output is redirected to  _log.txt_ and any error is redirected to the file _error.log_ inside the patient's folder.


```python
import os
import sys
import getopt
from subprocess import call
import re


def main():
    folder_name = '/home/mariano/DATA/Subtraction/'

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

<a name="ref1">[1]</a> M. Cabezas, A. Oliver, E. Roura, J. Freixenet, J.C Vilanova, Ll. Ramió-Torrentà, A. Rovira, X. Lladó. _Automatic multiple sclerosis lesion detection in brain MRI by FLAIR thresholding_. __Computer Methods and Programs in Biomedicine__, 115(3), pp. 147-161. 2014

<a name="ref2">[2]</a> M. Cabezas, J.F. Corral, A. Oliver, Y. Diez, M. Tintore, C. Auger, X. Montalban, X. Lladó, D. Pareto, A. Rovira. _Automatic multiple sclerosis lesion detection in brain MRI by FLAIR thresholding_. __American Journal of Neuroradiology__, to appear. 2016 
