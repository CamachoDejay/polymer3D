# README #

This repository holds all of the codes used in [insert Link to paper]. The way the data is processed is entirely described in the supplementary information of the articles.

### What is this repository for? ###

Tracking particles in three-dimension. The software allows for tracking particle in three-dimension when a time-series of z-stack is provided to it. If only a single 2D image is provided for each time point, the tracking is then perform in 2D.

### How do I get set up? ###

You will need to install bio-format for some of the extension (.his) see https://www.openmicroscopy.org/bio-formats/downloads/ where you can download the matlab toolbox.

#### 2D time-series ####
If you're data is composed of a time-series of 2D images (standard widefield or confocal microscope), the main code you will be using is the *mainTrackingExp.m*.

The first block two blocks of code in  *mainTrackingExp.m* aim are intended for user input:

path2ZCal and path2SRCal are only used for 3D analysis and need specialized calibration files, thus they should be left empty for the software to know that they should not be used:
path2ZCal = [];
path2SRCal = [];

file.path should be the path to the main folder of your analysis, that is, a folder that contains subfolder (e.g. mov1,mov2,mov3) which contains the data files. This convoluted system is due to the way our microscope encapsulated each dataset in a new subfolder which makes it easier for us to process. However, this allows to process multiple data set simultaneously. Typically you want to analyze data in group acquired in the same conditions. Here is an example of folder in which you want to process the folder containing data acquired on the 28 of August (the date, month, etc... are not needed, it is just an example taken out of the way I organize my data but there is no restriction as to the naming of the folders and subfolders as long as they exist on your computer or connected hard drive):

file.path  = 'C:\Users\X\Documents\2020-Data\August\28\'

file.ext should give the extension you are using. Currently **the code can only process .ome.tif (3D), .his and .tif files.** If you have another extension, we advice you to either convert the movie to .tif using imagej or to contact us so we can perform the modifications needed(see contact info at the bottom).

path2Cal is also for 3D multiplane data and should therefore be left empty here:
path2Cal = [];

The rest of the first block is linked to the detection and tracking parameters:
%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 80;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 500;
trackParam.memory  = 3;

We advice to keep the default parameter for most of these. The parameter that you might want to change is the chi2 which is sort of our detection threshold. If you are looking at very bright object going slightly above 80 could be good 80-120 (if you see that detection is done on diffraction rings for instance). On the other hand you might want to put it as low as 24 if you are working with very dim particles (e.g. single-molecules).

Regarding tracking the radius parameter describe the maximum distance you expect a particle have move between to consecutive time point. This value typically needs to be adapted as it will depend on both your framerate and the speed at which your object moves. Keep in mind that if a large value is picked there will still be some discrimination step if more than one particle is within the radius chosen while you can easily "lose" particles if you chose this radius too small, therefore we advice to always aim a bit higher than what you expect th motion to be.

The second block contains information about the processing of the data:

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.calibrate = false; %true to recalibrate;

info.type indicates wheter, in the movie, the object of interest are bright (darkfield or fluorescence == normal) or dark (transmission)
info.runMethod, indicates if you want to rerun the analysis(run) or load from potentially existing data (load). In the later case, the analysis will simply run if no data is found and thus is typically the default value. However, if you are not happy with the results you can rerun the analysis by changing this parameter to run.

info.frame2Load indicates the range of frame you want to process (e.g. 1:100, 23:75 or 'all')
info.fitMethod is what fit you want to use either the phasor approach from Martens et al 2018 or a 2D gaussian. The phasor approach is significantly faster and has similar accuracy so we picked it as default

info.zMethod is only used for 3D data so you do not need to mind its value in this case.
info.calibrate, false is default. Since the data is not 3D you most likely will not need to put it ever

Once these values are adapted to your need, you can press run ! A dialog should popup asking about the pixel size and the Full width half maximum you expect your particle of interest to have (in pixels)

The code will then run and hopefully after some minutes you will get a popup saying that your data was sucessfully saved. All the track traces from the different movies you analyzed will be stored in the main folder you provided in file.path in a single mat.file.

#### 3D time-series ####

We only use .ome.tif files and can only process multiplane microscopy data at the moment so if you are not using the same format we advice you to contact us (see contact below) for information about how to get setup or if you want us to extend the code so you can process your data. 

If using the same format as us here is the information you need to get setup (please read the 2D section too, only the info specific to 3D will be explained here):

In addition to the information provided above you will need to take a bit more care with the following:
info.zMethod is intensity by default and use the same methods as in the paper in which we used this code ([insert Link]). PSFE is possible but need an additional calibration file (zCalibration). 3DFit is still in development and is therefore not advice.

info.calibrate can be change to true if it turned out that your 2D calibration (that determines the distance between your planes, more details later) was not good and you recalculated, you will need to recalibrate the experimental data with the new calibration.

#### 2D Calibration ####
As explained in the paper we use a sample where fluorescent beads are spun cast onto a coverslip and we acquired a few z stack making sure that the range of the z stack is large enought so the beads will be in focus at some point in each of the planes. 

To get the 2D calibration file you will have to run *mainPlaneCalibration.m*

To run this code you just have to provide the main folder of analysis (folder containing subfolder each containing 1 Movie) for example:
file.path  = 'C:\Users\X\Documents\2020-Data\August\28\2DCal'
file.ext = '.ome.tif'

The data will be saved in the same path. When using *mainTrackingExp.m* you will need to provide that path as the calibraiton:
path2Cal = 'C:\Users\X\Documents\2020-Data\August\28\2DCal'

The software will automatically search for the calibration file that is saved by *mainPlaneCalibration.m*
 
#### SR Calibration ####

The SRCalibration is run on the same data as the 2D Calibration. You need to run *mainSRCalibration.m*. This code will correct mismatch between the different z planes with subpixel accuracy it is therefore not mandatory but adviced to increase the accuracy of your tracked trace.

To run this code you just have to provide the main folder of analysis (folder containing subfolder each containing 1 Movie) for example:
file.path  = 'C:\Users\X\Documents\2020-Data\August\28\2DCal'
file.ext = '.ome.tif'

The data will be saved in the same path. When using *mainTrackingExp.m* you will need to provide that path as the calibraiton:
path2SRCal = 'C:\Users\X\Documents\2020-Data\August\28\2DCal'

The software will automatically search for the SRcalibration file that is saved by *mainSRCalibration.m*.

Therefore, in *mainTrackingExp.m* both path2Cal and path2SRCal should have the same path given

#### Z Calibration ####
**This should only be used if using PSFE as a zMethod otherwise leave empty as described in the 2D**. To calculate a zCalibration, you will need to run mainZCalibration.m. You will need to have zStack of data suspended in 3D in a matrix as close as possible as the one you want to analyze but with the beads immobile. Potentially this could be run on the 2D calibration data but you risk to take into account abberation and optical effect from the 2D sample which might be different in your sample. The aim of this code is to calibrate the effect of the PSF engineering plate you have put in the system by quantifying the ellipticity as a function of z position.

We advice to left the parameter as default, you only need to provide the file.path and path2Cal since the data needs to be 2D calibrated to be thoroughly analyzed.

The path you give in File.path in *mainZCalibration* needs to be place in path2ZCal in *mainTrackingExp.m* after succesful Zcalibration.

Once you perform all the calibrations, you can follow the 2D part for the use of *mainTrackingExp.m* ! Remember that only the 2D Calibration is really mandatory and that is valid only for multiplane data. Confocal or other type of microscopy do not need to calibrate the position of the planes since you input it as a zstack parameter but we cant handle that type of data at the moment, please contact us for help.


### Contribution guidelines ###

///////////////////////////////

### Who do I talk to? ###

If you have any problem just contact Boris Louis (Boris.louis@kuleuven.be) or Rafael Camacho (rafael.camacho@gu.se)
