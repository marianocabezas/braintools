# BrainTools: C++ code based on ITK for brain MRI processing.

__BrainTools__ is a C++ repository with several fully automated tools to process brain MRI data from MS patients. Currently it implements _bias correction_ using the __N4__ algorithm; _atlas registration_ using affine and b-splines transformations; _tissue segmentation_ based on [1](https://github.com/NIC-VICOROB/braintools/edit/master/README.md#1); _lesion segmentation and detection_ also based on [1](1); and, _longitudinal lesion analysis_ based on [2].

Important:
This code supports version **3.20** of ITK. It has not been coded yet for newer versions and it won't compile with them. Also, the CMakelists might need some tweaking to define the ITK path depending on your OS.

[1] M. Cabezas, A. Oliver, E. Roura, J. Freixenet, J.C Vilanova, Ll. Ramió-Torrentà, A. Rovira, X. Lladó. _Automatic multiple sclerosis lesion detection in brain MRI by FLAIR thresholding_. __Computer Methods and Programs in Biomedicine__, 115(3), pp. 147-161. 2014

[2] M. Cabezas, J.F. Corral, A. Oliver, Y. Diez, M. Tintore, C. Auger, X. Montalban, X. Lladó, D. Pareto, A. Rovira. _Automatic multiple sclerosis lesion detection in brain MRI by FLAIR thresholding_. __American Journal of Neuroradiology__, to appear. 2016 
