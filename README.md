Matlab scripts for Direct Reconstruction in DCE-MRI
============================================================

(c) *Yi Guo*, *Sajan Goud Lingala*, *Krishna Nayak*, Jul 2016.

[Magnetic Resonance Engineering Laboratory](https://mrel.usc.edu)

**University of Southern California**

Code Structure
--------------
### Demo scripts
**Kt_Vp_SEN_AD_3d.m**
Reconstruction script to demonstrate direct reconstruction from under-sampled k-space in retrospective studies.
Please download data set at https://drive.google.com/file/d/0B4nLrDuviSiWajFDV1Frc3cxR0k/view?usp=sharing

### Functions: 
**conc2Ktrans.m**: 
Convert contrast concentration to TK maps
**conc2sig.m**: 
Convert contrast concentration to signal (images)
**genRGA.m**: 
Generate randomized golden-angle radial sampling pattern
**Kt_Vp_SEN.m**: 
Alternatively reconstruct Ktrans and Vp maps using l-bfgs
**Ktrans2conc.m**: 
Convert Ktrans maps to contrast concentration
**sig2conc2.m**: 
Convert signal (images) to contrast concentration
**Ktrans2sig_sen_WT.m**: 
cost and gradient calculation for the objective function for Ktrans
**Vp2sig_SEN_WT.m**: 
cost and gradient calculation for the objective funciton for Vp
**SAIF_p.m: 
generate population-averaged AIF
