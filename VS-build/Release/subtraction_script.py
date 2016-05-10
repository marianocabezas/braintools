import os
import sys
import getopt
from subprocess import call#, checkoutput
import re


def main():
    folder_name = 'C:\\Subtraccio\\'

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
        patient_folder = '%s%s\\' % (folder_name, patient)
        patient_done = '%sreport.html' % (patient_folder)
        if not os.path.isfile(patient_done):
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
                        timepoint_folder = '%s%s\\' % (patient_folder, timepoint)
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
                                    call([
                                        'ROBEX.exe',
                                        pd_image,
                                        pd_image_stripped,
                                        mask_image
                                    ], stdout=f, stderr=err)
                                    os.remove(pd_image_stripped)
                for baseline, followup in zip(timepoints[:-1], timepoints[1:]):
                    baseline_folder = '%s%s\\' % (patient_folder, baseline)
                    followup_folder = '%s%s\\' % (patient_folder, followup)
                    call([
                        'PreTool.exe',
                        baseline_folder,
                        followup_folder
                    ], stdout=f, stderr=err)
                    call([
                        'CoregTool.exe',
                        baseline_folder
                    ], stdout=f, stderr=err)
                    call([
                        'CoregTool.exe',
                        followup_folder
                    ], stdout=f, stderr=err)
                    call([
                        'AtlasTool.exe',
                        baseline_folder
                    ], stdout=f, stderr=err)
                    call([
                        'AtlasTool.exe',
                        followup_folder
                    ], stdout=f, stderr=err)
                    call([
                        'TissueTool.exe',
                        baseline_folder
                    ], stdout=f, stderr=err)
                    call([
                        'TissueTool.exe',
                        followup_folder
                    ], stdout=f, stderr=err)
                    call([
                        'SubtractionTool.exe',
                        baseline_folder,
                        followup_folder
                    ], stdout=f, stderr=err)
                    call([
                        'DeformableTool.exe',
                        baseline_folder,
                        followup_folder
                    ], stdout=f, stderr=err)
                    call([
                        'PostTool.exe',
                        baseline_folder,
                        followup_folder
                    ], stdout=f, stderr=err)
                with open(patient_done, 'w') as report:
                    pass
                    # HTML header and other stuff
                    # Volumes and number of lesions
                    #for followup in timepoints[1:]:
                        #subtraction_folder = '%s%s\\' % (patient_folder, followup)
                        #report =
                        #output = check_output([
                        #    'AnalyseTool.exe',
                        #    subtraction_folder
                        #], stderr=err)
                        #parsed_output = re.search(r'(\w|\W)*;(\w|\W)*;(\w|\W)*;',file)
                        #parsed_output.mat

if __name__ == '__main__':
    main()