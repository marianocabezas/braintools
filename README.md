# BrainTools: C++ code based on ITK for brain MRI processing.

__BrainTools__ is a C++ repository with several fully automated tools to process brain MRI data from MS patients. Currently it implements _bias correction_ using the __N4__ algorithm; _atlas registration_ using affine and b-splines transformations; _tissue segmentation_ based on [1]; _lesion segmentation and detection_ also based on [1]; and, _longitudinal lesion analysis_ based on [2].

Important:
This code supports version **3.20** of ITK. It has not been coded yet for newer versions and it won't compile with them. Also, the CMakelists might need some tweaking to define the ITK path depending on your OS.


