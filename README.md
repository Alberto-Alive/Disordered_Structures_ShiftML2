# Scripts overview

  ## Make predictions
  **Predict the chemical shift for a set of xyz files containing your chemical structures**
  - outputs a magres file for each xyz file and a csv file containing all the predictions run on the xyz files
  
  ## Output spectra
  **Create a spectrum for each magres file**
  - outputs a png file for each magres file
  
  ## Compare spectra
  **Calculate the non/weighted average spectrum across all magres files and identify similar spectra using a heatmap and a plot of subplots**
  - outputs a png file for the non/weighted average spectrum and for the plot of stacked spectra  and for the heatmap
  
 # Setup
 
  ## 1. Python environment
  
  Windows
  Install WSL2 (Windows Subsystem for Linux 2). The scripts were tested using Ubuntu 22.042.2 LTS.
  Create a conda environment where you install all the dependencies (ShiftML2, Soprano, ASE, Click, Matplotlib, Numpy, Pandas, Seaborn etc).
      
  ## 2. Get ShiftML2 models
  ShiftML2 models are pretrained models for predicting chemical shifts for specific chemical elements: C, Ca, Cl, F, H, K, Mg, N, Na, O, P, S.
  Download them from: [Zenodo](https://zenodo.org/record/6782654) or use this doi to find it: 10.5281/zenodo.6782654
   
 # Run scripts
  
  ## 1. Make predictions: make_predictions.py
  
    python make_predictions.py *.xyz ./example/ababub_xyz_files/*.xyz H -ase "{'index' : ':'}" -s
  
  ## 2. Output spectra: output_spectra.py
    python output_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -fr "{'min':10, 'max':35}"
  
  ## 3. Compare spectra: compare_spectra.py
    python compare_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -w [1, 4, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1] -fr "{'min':10, 'max':35}"
