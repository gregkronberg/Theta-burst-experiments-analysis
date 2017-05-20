Synaptic Plasticity Analysis

1. Original data are exported from LabChart as .mat files and stored in a 
	"Raw Data" directory.
	Images are saved in the "Raw Images" directory

2. tbs_preprocessing.m files search the Raw Data directory for new files.
	If new files are found, the raw data are organized into matrices 
	and saved in the Processed Data directory.  Information about 
	the experimental parameters is extracted from comments made in LabChart
	and stored in a structure called slices.  slices is a structure that 
	contains information about all processed slices, organized as
	slices{conditions}(slice_number).parameter
	slices is saved with a "conditions" variable that lists all conditions.
	Images are preprocessed separately, using the tbs_preprocessing_images.m 
	file

3. tbs_vargen.m files use the processed data to create global variables for further analysis.
	For example fEPSP area during induction for all slices is stored as
	a variable structure called fepsp_area in the Matlab Variables folder.
	All variables are organized like the slices structure as 
	variable{conditions}(slice_number).data, where slice_number corresponds
	to slice number in slices

4. tbs_analysis.m files then perform further analyses after the variable generation stage.
	The principle difference between the analysis and vargen scripts is that vargen scripts
	requiring looping over and loading all processed data files, which include large time
	series files.  Analysis files can run quicker by loading global variables and running further
	analyses.  Analysis variables are added to the global variable structures that they were derived 
	from and the global variables are saved 
	


Image Data