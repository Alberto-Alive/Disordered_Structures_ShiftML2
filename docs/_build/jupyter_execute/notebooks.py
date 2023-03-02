#!/usr/bin/env python
# coding: utf-8

# # Scripts
# 
# The package contains 3 scripts for: making predictions, creating spectra, and comparing spectra. 
# 

# ## Make predictions
# Predict the chemical shift for a set of xyz files containing your chemical structures.
# 
# Outputs a magres file for each xyz file and a csv file containing all the predictions run on the xyz files.

# ```
# import os
# from datetime import datetime
# 
# import ase
# import ase.io
# import numpy as np
# import click
# 
# import ShiftML2 as sml
# import time
# ```

# ```
# def validate_input(ase_params, xyz_files_paths, chemical_element):
#     for path in xyz_files_paths:
#         if not path.endswith(".xyz"):
#             raise click.BadParameter(f"{path} is not a valid path targetting an .xyz file.")
#     if ase_params:
#         try:
#             ase_params = eval(ase_params)
#         except:
#             raise TypeError(f"You must provide the dict as a string using double quotes\nTry using: \"{ase_params}\"") 
#         if not isinstance(ase_params, dict):
#             raise TypeError("The ase params is not a dictionary\nProvide an in double quotes dictionary of the form: \"{'index' : ':'}\"")
#     
#     if not isinstance(chemical_element, str):
#         raise TypeError(f"The following {chemical_element} is not a string. Please provide just one chemical element as string.")
#     elif chemical_element not in ["H", "C", "N", "O", "Cl", "Ca", "F", "K", "Mg", "N", "Na", "P", "S"]:
#         raise KeyError(f"The following chemical element '{chemical_element}' is not a supported chemical element:[H, C, N, O, Cl, Ca, F, K, Mg, N, Na, P, S].")
#     return ase_params 
# ```

# ```
# @click.command()
# @click.argument("xyz_files_paths", nargs=-1, type=click.Path(exists=True))
# @click.argument('chemical-element')
# @click.option("--ase_params", '-ase', default={'index' : ':'}, help="ASE parameters for the structure. The dictionary should have the form: \"{'index' : ':'}\"")
# @click.option("--start", '-s', is_flag=True, help="Start the script.")
# 
# def cli(chemical_element, xyz_files_paths, ase_params, start):
#     ase_params = validate_input(ase_params = ase_params, xyz_files_paths = xyz_files_paths, chemical_element = chemical_element)
#     structure = ChemicalStructure(chemical_element, xyz_files_paths, ase_params)
#     if start:
#         result = structure.start_script()
#         print(f"Result of prediction: {result}")
#         print("Done!")
#     else:
#         print(f"Following structure {structure} was initialised. Use -s flag to start the script.")
# ```
# 

# ```
# class ChemicalStructure(object):
#     """ChemicalStructure class
# 
#        A class representing the structure of a molecule.This class will
#        run ShiftML2 on all the xyz files found in the provided paths and 
#        predict the chemical shifts of all configurations in the file 
#        (unless user specified using the ase_params).
# 
#        | Args:
#        |   xyz_files_paths (list): list of xyz file paths containing the chemical
#        |                     structures the chemical shift prediction is made for.
#        |   chemical_element (str): letter symbol of the chemical element to have
#        |                      the chemical shift predictions applied to: C, Ca, Cl,
#        |                      F, H, K, Mg, N, Na, O, P, S.
#        |   ase_params (dict): ASE read function parameters (ase.io.read), default
#        |                   is ':', returning all configurations from the structure.
#        |
#     Example:
#     >>> python ChemicalStructure.py *.xyz ./fake_folder/*.xyz H -ase "{'index' : ':'}" -s
#     """
# 
# 
# 
#     def __init__(self, chemical_element, xyz_files_paths=None, ase_params={'index' : ':'}):
#         self.xyz_files_paths = xyz_files_paths
#         self.magres_array_lengths = []
#         self.arrays_for_magres = []
#         self.structs_list = []
#         self.names_list = []
#         self.ase_params = ase_params
#         self.chemical_element = chemical_element
#         self.structs = None
#         self.y_pred =None
#         self.y_err = None
#         self.model = None
#         self.csv_file_name = f"{chemical_element}-{(datetime.now().isoformat(sep='_', timespec='seconds'))}-ShiftML2-prediction"
#     
# 
#     def start_script(self):
#         """Start script function 
#         
#             A function that makes the calls to the appropriate class functions
#             to make the chemical shift predictions and output the csv and magers files.
# 
#         """
# 
#         self.read_xyz_files(self.xyz_files_paths)
#         self.load_model()
#         self.build_model_representation()
#         self.predict_set()
#         self.output_csv()
#         self.prepare_magres()
#         self.output_magres()
# 
# 
#     def read_xyz_files(self, list_of_paths):
#         """Read xyz file path
# 
#             A function that reads all the xyz files in an array of xyz file paths.
#             The function creates a list of atom objects collected from all xyz files.
#             It also creates an list of file names 'names_list' used to name the magres files.
# 
#         | Args:
#         |   Input:
#         |       list_of_paths (list): list of paths leading to xyz files.
#         |
#         |   Output:
#         |       structs_list (list): list of atom objects found in every xyz file.
#         |       names_list (list): List of xyz files of format: 
#         |       {xyz_file_name}-{file_index}-atom_obj-{atom_object_index}-ShiftML2-prediction.
#         |       magres_array_lengths (list): List of lengths corresponding to
#         |                               number of time the chemical element 
#         |                               to be predicted, appears in the atom 
#         |                               objects found within an xyz file.
#         |         
#         
#         """
#         for f_index, xyz_file_path in enumerate(list_of_paths):
#             print("Path: ",xyz_file_path)
#             try:
#                 structs = ase.io.read(xyz_file_path, **self.ase_params)
#             except:
#                 raise ValueError(f"Failed to read file. \nWrong 'ase_params' provided: {self.ase_params}")
#             if not structs:
#                 raise ValueError("No structures could be read from the paths provided.")
#             file_name = os.path.basename(xyz_file_path)
#             index_of_dot = file_name.index('.')
#             file_name_without_extension = file_name[:index_of_dot]
#             if type(structs) == list:
#                 for o_index, atom_object in enumerate(structs):
#                     self.structs_list.append(atom_object)
#                     self.names_list.append(f"{file_name_without_extension}-{f_index:04d}-atom_obj-{o_index:04d}")
#                     target_atoms = [atom for atom in atom_object if atom.symbol == self.chemical_element]
#                     num_atoms = len(target_atoms)
#                     self.magres_array_lengths.append(num_atoms) 
#             elif type(structs) == 'ase.atoms.Atoms':
#                 self.structs_list.append(structs)
#                 self.names_list.append(f"{file_name}-{f_index:04d}-atom_obj-0000")
#                 target_atoms = [atom for atom in structs if atom.symbol == self.chemical_element]
#                 num_atoms = len(target_atoms) 
#                 self.magres_array_lengths.append(num_atoms)  
#             else:
#                 raise TypeError(f"The type of {structs} should be either 'list' or 'ase.atoms.Atoms'")
# 
#  
#     def load_model(self):
#         """Get model function
# 
#             Get ShiftML2 model for specific chemical element.
# 
#         | Args:
#         |   Input:
#         |       chemical_element (str): Element to predict ("H", "C", "N", "O", "S", 
#         |                                 "F", "P", "Cl", "Na", "Ca", "Mg" or "K").
#         |                    
#         |   Output:
#         |       model (ShiftML2.ShiftML2.Model): ShiftML2 model class.
#         | 
# 
#         """
#         self.model = sml.Model(self.chemical_element)
# 
# 
#     def build_model_representation(self):
#         """Build model function
# 
#             Build the descriptor representation of atomic environments of a given atom in a set of structures.
# 
#         | Args:
#         |   Input:
#         |       structs_list (list): List of atom objects found in every xyz file path.  
#         |                    
#         |   Output:
#         |       model (ShiftML2.ShiftML2.Model): ShiftML2 model class.
#         | 
# 
#         """
#         self.model.build_representation(self.structs_list)
# 
# 
#     def predict_set(self):
#         """Real prediction function
# 
#             Predict target values on a test set using an ensemble of trained models.
# 
#         | Args:
#         |   Input:
#         |       structs_list (list): List of atom objects found in every xyz file path.  
#         |                    
#         |   Output:
#         |       model (ShiftML2.ShiftML2.Model): ShiftML2 model class.
#         |       y_pred2 (numpy.ndarray): Predicted values for the real set.
#         |       y_err2 (numpy.ndarray): Predicted errors for the real set.
#         | 
# 
#         """
#         start_time = time.time()
#         self.y_pred, self.y_err = self.model.predict()
#         end_time = time.time()
#         elapsed_time = end_time - start_time
#         print("Elapsed time: ", elapsed_time)
# 
#             
#     def output_csv(self):
#         """Output csv function
# 
#             Output one csv file for all xyz files found in user provided paths.
# 
#         | Args:
#         |   Input:
#         |       y_pred2 (numpy.ndarray): Predicted values for the real set.
#         |       y_err2 (numpy.ndarray): Predicted errors for the real set.
#         |
#         |   Output:
#         |       file (csv): the csv file will be stored within the current directory.
#         |                                   
#         | 
# 
#         """
#         prep_array = np.asarray([self.y_pred, self.y_err])
#         np.savetxt(f"./{self.csv_file_name}.csv", prep_array, delimiter=",")
# 
# 
#     def prepare_magres(self):
#         """Prepare magres function
# 
#             Slice predicted values returned as one big array into smaller arrays to be
#             replaced their corresponding values at corresponding positions within their 
#             respective atom object and transform it into a magres file.
# 
#         | Args:
#         |   Input:
#         |      y_pred2 (numpy.ndarray): Predicted values for the real set.  
#         |      magres_array_lengths (list): List of lengths corresponding to
#         |                               number of time the chemical element 
#         |                               to be predicted, appears in the atom 
#         |                               objects found within an xyz file.
#         |
#         |   Output:
#         |       array_for_magres (list): List of predicted values for the real set 
#         |          corresponding to their atom object found in a specific xyz file.     
#         |       arrays_for_magres (list): List of array_for_magres arrays.
#         | 
# 
#         """
#         start = 0
#         for slice_length in self.magres_array_lengths:
#             end = start + slice_length
#             array_for_magres = self.y_pred[start:end]
#             self.arrays_for_magres.append(array_for_magres)
#             start = end
# 
# 
#     def output_magres(self):
#         """Output magres function
# 
#             Output one magres file for each atom obect found in every xyz file.
# 
#         | Args:
#         |   Input:
#         |       arrays_for_magres (numpy.ndarray): List of lists of predicted values
#         |                                          for the real set.  
#         |       structs_list (list): List of atom objects found in every xyz file path.
#         |       names_list (list): List of xyz files of format: 
#         |       {xyz_file_name}-{file_index}-atom_obj-{atom_object_index}-ShiftML2-prediction.
#         |       chemical_element (str): Element to predict ("H", "C", "N", "O", "S", 
#         |                                 "F", "P", "Cl", "Na", "Ca", "Mg" or "K").
#         |       indices (list): label column found in xyz. if label column if missing
#         |                 it will be replaced with their corresponding chemical element symbol.
#         |
#         |   Output:
#         |       file (magres): magres files will be stored within the current directory.
#         |                                   
#         | 
# 
#         """
#         for i_list, (pred_array, atoms, name) in enumerate(zip(self.arrays_for_magres, self.structs_list, self.names_list)):
#             ms = np.zeros((len(atoms),3,3))
#             atom_index = 0
#             for i_atom, atom in enumerate(atoms):
#                 if atom.symbol == self.chemical_element:
#                     ms[i_atom][np.diag_indices(3)] = pred_array[atom_index]
#                     atom_index += 1
#             atoms.arrays['ms'] = ms
#             if atoms.has('labels'):
#                 labels = list(atoms.get_array('labels'))
#             else:
#                 labels = atoms.get_chemical_symbols()  
#             indices = [labels[:i + 1].count(labels[i]) for i in range(len(labels))]
#             atoms.arrays['indices'] = indices
#             atoms.info['magres_units'] = {'ms': 'ppm'}
#             atoms.write(f'./{name}-ShiftML2-prediction.magres', format='magres')
#         print(f"Outputted magres files count: {i_list + 1}")
# 
# if __name__ == '__main__':
#     cli()
# ```

# ## Create spectra
# Create a spectrum for each magres file.
# 
# Outputs a png file for each magres file.
# 

# ```
# import os
# 
# import numpy as np
# import click
# 
# from ase.io import read
# from soprano.calculate.nmr import NMRCalculator
# from soprano.properties.nmr import MSIsotropy
# import matplotlib.pyplot as plt
# ```

# ```
# def validate_input(title, chemical_element, path_to_magres, frequency_range ):
#     if not isinstance(title, str):
#         raise TypeError(f"The following {title} is not a string.")
#     elif len(title) > 30:
#         raise ValueError(f"The title provided is too long. Max characters permitted is 30")
#     
#     if not isinstance(chemical_element, str):
#         raise TypeError(f"The following {chemical_element} is not a string. Please provide just one chemical element as string.")
#     elif chemical_element not in ["H", "C", "N", "O", "Cl", "Ca", "F", "K", "Mg", "N", "Na", "P", "S"]:
#         raise KeyError(f"The following chemical element '{chemical_element}' is not a supported chemical element:[H, C, N, O, Cl, Ca, F, K, Mg, N, Na, P, S].")
#     
#     if frequency_range:
#         try:
#             frequency_range = eval(frequency_range)
#         except:
#             raise TypeError(f"You must provide the dict as a string using double quotes\nTry using: \"{frequency_range}\"") 
#         if not isinstance(frequency_range, dict):
#             raise TypeError("The frequency range is not a dictionary\nProvide an in double quotes dictionary of the form: \"{'min': 0, 'max': 1000}\"")
#     for index, path in enumerate(path_to_magres):
#         if not path.endswith(".magres"):
#             raise click.BadParameter(f"{path} is not a valid path targetting .magres files.")
#         else:
#             print(f"Path for plot {index:02d}: {path}")    
#     atoms_list = [read(path) for path in path_to_magres]
#  
#     print(f'Magres files count: {len(path_to_magres)}\n')
#     return atoms_list, frequency_range
# ```

# ```
# @click.command()
# @click.argument('path_to_magres', nargs=-1, type=click.Path(exists=True))
# @click.argument('chemical-element')
# @click.option('--freq_broad', '-fb', type=float, required=True, help='The frequency broadening factor to use. [required]')
# @click.option('--bins_count', '-b', type=int, default=300, help='The number of bins to use in the spectrum. If not provided, the default is 300.')
# @click.option('--reference', '-r', type=str, default=None, help='The reference atom for the chemical shifts.')
# @click.option('--title', '-t', type=str, default=None, help='A title for the plot.')
# @click.option('--frequency_range', '-fr', default=None, help='A dictionary specifying the minimum and maximum frequency to use in the spectrum. The keys should be "min" and "max".')
# ```

# ```
# def create_spectra(chemical_element,  freq_broad, bins_count, path_to_magres=None, reference=None, title=None, frequency_range=None):
#     """
#     This function generates a 1D powder NMR spectrum for the given chemical element in the `.magres` files found in the specified path.
#     
#     Parameters:
#     path_to_magres (str, optional): The path to the `.magres` file(s) to use for the calculation.
#     chemical_element (str): The chemical element to generate the spectrum for (e.g. "H").
#     freq_broad (float): The frequency broadening factor to use.
#     bins_count (int): The number of bins to use in the spectrum.
#     reference (str, optional): The reference atom for the chemical shifts.
#     title (str, optional): A title for the plot.
#     frequency_range (dict, optional): A dictionary specifying the minimum and maximum frequency to use in the spectrum. The keys should be "min" and "max".
# 
#     Returns:
#     None
# 
#     Raises:
#     ValueError: If the specified path does not lead to a `.magres` file or a folder containing `.magres` files.
# 
#     Example:
#     >>> python create_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -fr "{'min':0, 'max':1000}"
#     """
#     atoms_list, frequency_range = validate_input(title = title, chemical_element = chemical_element, path_to_magres = path_to_magres, frequency_range = frequency_range )
#     print('Plotting the plots..')
#     max_freq = None
#     min_freq = None
#     for index, (atoms, path) in enumerate(zip(atoms_list, path_to_magres)):
#         print(f"Start plotting: plot {index:02d}")
#         calc = NMRCalculator(atoms)
# 
#         if reference:
#             calc.set_reference(reference, chemical_element)
#             use_reference = True
#         else:
#             use_reference = False
# 
#         calc.set_powder(N=1)
# 
#         if not frequency_range:
#             iso = MSIsotropy.get(atoms)
#             max_freq=max(iso)*1.1
#             min_freq=min(iso)*1.1
#             spec_elem_iso, freqs = calc.spectrum_1d(chemical_element, max_freq=max_freq, min_freq=min_freq, bins=bins_count,
#                                             freq_broad=freq_broad, use_reference=use_reference)
#         else:
#             spec_elem_iso, freqs = calc.spectrum_1d(chemical_element, max_freq=frequency_range['max'], min_freq=frequency_range['min'], bins=bins_count,
#                                             freq_broad=freq_broad, use_reference=use_reference)
#         spec_elem_iso = spec_elem_iso / np.sum(spec_elem_iso)
# 
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
# 
#         ax.plot(freqs, spec_elem_iso, label='Isotropic shifts')
#         if reference:
#             ax.set_xlabel(f'{chemical_element} chemical shift (ppm)')
#             ax.invert_xaxis()
#         elif reference == None:
#             ax.set_xlabel(f'{chemical_element} chemical shielding (ppm)')
#         ax.set_ylabel('Intensity (a.u.)')
#         if title:
#             ax.set_title(title)
# 
#         if frequency_range and not reference:
#             ax.set_xlim(frequency_range['min'], frequency_range['max'])
#         elif frequency_range and reference:
#             ax.set_xlim(-frequency_range['min'], -frequency_range['max'])
#         elif not frequency_range and not reference: 
#             ax.set_xlim(min_freq, max_freq)
#         elif not frequency_range and reference:
#             ax.set_xlim(-min_freq, -max_freq)
#         ax.set_xlim()
#         ax.set_yticks([])
#         ax.legend()
# 
#         file_name = os.path.basename(path)
#         index_of_dot = file_name.index('.')
#         file_name_without_extension = file_name[:index_of_dot]
#         plt.savefig(f'./{file_name_without_extension}-spectrum.png', dpi=300)
#         plt.close()
#         print(' Finished plotting.')
#     print("Done!")
# if __name__ == '__main__':
#     create_spectra()
# ```

# ## Compare spectra
# Calculate the non/weighted average spectrum across all magres files and identify similar spectra using a heatmap and a plot of subplots.
# 
# Outputs a png file for the non/weighted average spectrum and for the plot of stacked spectra and for the heatmap.

# ```
# import numpy as np
# import seaborn as sns
# 
# from ase.io import read
# from soprano.properties.nmr import MSIsotropy
# from soprano.calculate.nmr import NMRCalculator
# 
# import matplotlib.pyplot as plt
# import click
# ```

# ```
# def validate_input(title, chemical_element, path_to_magres, frequency_range, weights):
#     if not isinstance(title, str):
#         raise TypeError(f"The following {title} is not a string.")
#     elif len(title) > 30:
#         raise ValueError(f"The title provided is too long. Max characters permitted is 30")
#     
#     if not isinstance(chemical_element, str):
#         raise TypeError(f"The following {chemical_element} is not a string. Please provide just one chemical element as string.")
#     elif chemical_element not in ["H", "C", "N", "O", "Cl", "Ca", "F", "K", "Mg", "N", "Na", "P", "S"]:
#         raise KeyError(f"The following chemical element '{chemical_element}' is not a supported chemical element:[H, C, N, O, Cl, Ca, F, K, Mg, N, Na, P, S].")
# 
#     if frequency_range:
#         try:
#             frequency_range = eval(frequency_range)
#         except:
#             raise TypeError(f"You must provide the dict as a string using double quotes\nTry using: \"{frequency_range}\"") 
#         if not isinstance(frequency_range, dict):
#             raise TypeError("The frequency range is not a dictionary\nProvide an in double quotes dictionary of the form: \"{'min': 0, 'max': 1000}\"")
# 
#     if weights:
#         try:
#             weights = eval(weights)
#         except:
#             raise TypeError(f"You must provide the list as a string using double quotes\nTry using: \"{weights}\"")
#         if not isinstance(weights, list):
#             raise TypeError(f"Provided {weights} is not a list, use format: -w [0,1,1]")
#         elif weights and (len(path_to_magres) != len(weights) or min(weights) < 0):
#             raise ValueError(f"You need an arrray of weights (positive values including 0) that is of length {len(path_to_magres)} \nYour array of weights has a length of {len(weights)}:{weights}") 
#     for index, path in enumerate(path_to_magres):
#         if not path.endswith(".magres"):
#             raise click.BadParameter(f"{path} is not a valid path targetting .magres files.")
#         else:
#             print(f"Path for plot {index:02d}: {path}")    
#     atoms_list = [read(path) for path in path_to_magres]
#     
#     print(f'Magres files count: {len(path_to_magres)}\n')
#     return atoms_list, frequency_range, weights
# ```

# ```
# @click.command()
# @click.argument('path_to_magres', nargs=-1, type=click.Path(exists=True))
# @click.argument('chemical-element')
# @click.option('--title', '-t', type=str, required=True, help='Title of the average spectra. [required]')
# @click.option('--freq-broad', '-fb', type=float, default=0.05, help='Broadening applied to the spectra. If not provided, the default is 0.05.')
# @click.option('--bins-count', '-b', type=int, default=300, help='Number of bins in the spectra. If not provided, the default is 300.')
# @click.option('--weights', '-w', default=None, help='List of weights (N=[0, ∞]) for each spectra. If not provided, a mean is taken. Weights are not normalised, for example given the weights [0, 1, 1] spectra of index 0 will not count but spectra of index 1 and 2 will count 50% each.')
# @click.option('--reference', '-r',  type=float, default=None, help='Reference to use for shifting the spectra. If not provided, no reference is used.')
# @click.option('--frequency-range', '-fr', default=None, help='Dictionary specifying the maximum and minimum frequency range to use (use double quotes around the dict). If not provided, range is calculated using chemical element.')
# ```

# ```
# def compare_spectra(chemical_element, title, freq_broad, bins_count, weights=None, path_to_magres=None, reference=None, frequency_range=None):
#     """
#     This function computes the average spectra of chemical element over multiple ".magres" files.
#     
#     Parameters:
#     path_to_magres (str, required): Path to the folder containing the ".magres" files. Defaults to None.
#     chemical_element (str, required): Chemical element to analyze.
#     title (str): Title of the average spectra.
#     freq_broad (float): Broadening applied to the spectra.
#     bins_count (int): Number of bins in the spectra.
#     weights (list, optional): List of weights for each spectra. If not provided, a mean is taken. Defaults to None. Weights are not normalised, for example given the weights [0, 1, 1] spectra of index 0 will not count but spectra of index 1 and 2 will count 50% each.
#     reference (float, optional): Reference to use for shifting the spectra. If not provided, no reference is used. Defaults to None.
#     frequency_range (dict, optional): Dictionary specifying the maximum and minimum frequency range to use. If not provided, range is calculated using chemical element. Defaults to None.
#     
#     Returns:
#     None: The function saves the subplots of individual spectra and the average spectra as PNG files in the current directory.
#     
#     Raises:
#     ValueError: If the path to the ".magres" files is not provided, or if the path leads to less than 2 or more than 30 ".magres" files.
#     ValueError: If the length of weights array is not equal to the number of spectra, or if there are negative values in the weights array.
#     
#     Example:
#     >>> python compare_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -w [0,1,1] -fr "{'min':0, 'max':1000}"
#     >>> compare_spectra("title", "C", 0.05, 300, reference="N", frequency_range={"min":0, "max":1000})
#     """
#     atoms_list, frequency_range, weights = validate_input(title = title, chemical_element = chemical_element, path_to_magres = path_to_magres, frequency_range = frequency_range, weights = weights)
#     
#     #Plotting subplots
#     print('Plotting the subplots..')
#     iso_spectras = []
#     freq_spectras = []
#     max_freq = None
#     min_freq = None
#     for atoms in atoms_list:
#         calc = NMRCalculator(atoms)
# 
#         if reference:
#             calc.set_reference(reference, chemical_element)
#             use_reference = True
#         else:
#             use_reference = False
# 
#         calc.set_powder(N=1)
# 
#         if not frequency_range:
#             iso = MSIsotropy.get(atoms)
#             max_freq=max(iso)*1.1
#             min_freq=min(iso)*1.1
#             spec_elem_iso, freqs = calc.spectrum_1d(chemical_element, max_freq=max_freq, min_freq=min_freq, bins=bins_count,
#                                             freq_broad=freq_broad, use_reference=use_reference)
#         else:
#             spec_elem_iso, freqs = calc.spectrum_1d(chemical_element, max_freq=frequency_range['max'], min_freq=frequency_range['min'], bins=bins_count,
#                                             freq_broad=freq_broad, use_reference=use_reference)
#         spec_elem_iso = spec_elem_iso / np.sum(spec_elem_iso)
#         iso_spectras.append(spec_elem_iso)
#         freq_spectras.append(freqs)
# 
#     number_of_subplots = len(iso_spectras)
#     fig, axs = plt.subplots(nrows=number_of_subplots, sharex=True, figsize=(12, 22))
#     plt.subplots_adjust(hspace=1.5)
#     fig.suptitle("Subplots of magres", fontsize=18, y=0.95)
# 
#     for elem_freq, elem_iso, ax in zip(freq_spectras, iso_spectras, axs.ravel()):
#         ax.plot(elem_freq, elem_iso)
#         ax.set_ylabel('Intensity (a.u.)')
#     
#     if reference:
#         axs[-1].set_xlabel(f'{chemical_element} chemical shift (ppm)')
#         ax.invert_xaxis()
#     elif reference == None:
#         axs[-1].set_xlabel(f'{chemical_element} chemical shield (ppm)')
#     if frequency_range and not reference:
#         ax.set_xlim(frequency_range['min'], frequency_range['max'])
#     elif frequency_range and reference:
#         ax.set_xlim(-frequency_range['min'], -frequency_range['max'])
#     elif not frequency_range and not reference: 
#         ax.set_xlim(min_freq, max_freq)
#     elif not frequency_range and reference:
#         ax.set_xlim(-min_freq, -max_freq)
#     plt.savefig(f'./{title}-{number_of_subplots}-spectrum_subplots.png', dpi=300)
#     plt.close()
#     print(' Finished plotting.')
# 
#     #Averaging plot
#     print('Averaging subplots..')
#     if weights:
#         average_iso = np.average(iso_spectras, weights=weights, axis=0)
#     else:
#         average_iso = np.mean(iso_spectras, axis=0)
#     fig, ax = plt.subplots(figsize=(12, 8))
#     ax.plot(freqs, average_iso)
# 
#     if frequency_range and not reference:
#         ax.set_xlim(frequency_range['min'], frequency_range['max'])
#     elif frequency_range and reference:
#         ax.set_xlim(-frequency_range['min'], -frequency_range['max'])
#     elif not frequency_range and not reference: 
#         ax.set_xlim(min_freq, max_freq)
#     elif not frequency_range and reference:
#         ax.set_xlim(-min_freq, -max_freq)
# 
#     if reference:
#         ax.set_xlabel(f'{chemical_element} chemical shift (ppm)')
#         ax.invert_xaxis()
#     elif reference == None:
#         ax.set_xlabel(f'{chemical_element} chemical shielding (ppm)')
#     
#     ax.set_ylabel('Intensity (a.u.)')
#     fig.suptitle("Average of magres subplots", fontsize=18, y=0.95)
#     plt.savefig(f'./{title}-average-spectrum.png', dpi=300)
#     plt.close()
#     print(' Finished averaging.')
# 
#     #Correlation Matrix
#     print('Calculating correlation..')
#     corr_matrix = np.corrcoef(iso_spectras)
#     n = corr_matrix.shape[0]
#     groups = []
#     for i in range(n):
#         for j in range(i+1, n):
#             if corr_matrix[i, j] >= 0.99:
#                 groups.append((i, j))
# 
#     unmatched = list(range(n))
#     all_groups = []
#     while len(unmatched) > 0:
#         i = unmatched.pop(0)
#         matches = [i]
#         inds = [j for j in unmatched if corr_matrix[i][j] >= 0.99]
#         matches += inds
#         unmatched = [u for u in unmatched if u not in matches]
#         all_groups.append(matches)
# 
#     plot_title = ''
#     for group in all_groups:
#         plot_title += ' ' + str(group)
# 
#     plt.figure(figsize=(8, 6))
#     sns.heatmap(corr_matrix, annot=True, cmap="YlGnBu", annot_kws={"size": 8})
#     plt.xlabel("Spectrum index")
#     plt.ylabel("Spectrum index")
#     plt.title(plot_title)
#     plt.savefig(f'./{title}-correlation-matrix.png', dpi=300)
#     plt.close()
#     for i, fname in enumerate(path_to_magres):
#         print(i, fname)
#     print(' Finished calculation.')
# 
# if __name__ == '__main__':
#     compare_spectra()
# ```
