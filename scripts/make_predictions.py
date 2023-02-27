import os
from datetime import datetime

import ase
import ase.io
import matplotlib.pyplot as plt
import numpy as np
import click

import ShiftML2 as sml
import time

def validate_input(ase_params, xyz_files_paths, chemical_element):
    for path in xyz_files_paths:
        if not path.endswith(".xyz"):
            raise click.BadParameter(f"{path} is not a valid path targetting an .xyz file.")
    if ase_params:
        try:
            ase_params = eval(ase_params)
        except:
            raise TypeError(f"You must provide the dict as a string using double quotes\nTry using: \"{ase_params}\"") 
        if not isinstance(ase_params, dict):
            raise TypeError("The ase params is not a dictionary\nProvide an in double quotes dictionary of the form: \"{'index' : ':'}\"")
    
    if not isinstance(chemical_element, str):
        raise TypeError(f"The following {chemical_element} is not a string. Please provide just one chemical element as string.")
    elif chemical_element not in ["H", "C", "N", "O", "Cl", "Ca", "F", "K", "Mg", "N", "Na", "P", "S"]:
        raise KeyError(f"The following chemical element '{chemical_element}' is not a supported chemical element:[H, C, N, O, Cl, Ca, F, K, Mg, N, Na, P, S].")
    return ase_params 

@click.command()
@click.argument("xyz_files_paths", nargs=-1, type=click.Path(exists=True))
@click.argument('chemical-element')
@click.option("--ase_params", '-ase', default={'index' : ':'}, help="ASE parameters for the structure. The dictionary should have the form: \"{'index' : ':'}\"")
@click.option("--start", '-s', is_flag=True, help="Start the script.")

def cli(chemical_element, xyz_files_paths, ase_params, start):
    ase_params = validate_input(ase_params = ase_params, xyz_files_paths = xyz_files_paths, chemical_element = chemical_element)
    structure = ChemicalStructure(chemical_element, xyz_files_paths, ase_params)
    if start:
        result = structure.start_script()
        print(f"Result of prediction: {result}")
        print("Done!")
    else:
        print(f"Following structure {structure} was initialised. Use -s flag to start the script.")

class ChemicalStructure(object):
    """ChemicalStructure class

       A class representing the structure of a molecule.This class will
       run ShiftML2 on all the xyz files found in the provided paths and 
       predict the chemical shifts of all configurations in the file 
       (unless user specified using the ase_params).

       | Args:
       |   xyz_files_paths (list): list of xyz file paths containing the chemical
       |                     structures the chemical shift prediction is made for.
       |   chemical_element (str): letter symbol of the chemical element to have
       |                      the chemical shift predictions applied to: C, Ca, Cl,
       |                      F, H, K, Mg, N, Na, O, P, S.
       |   ase_params (dict): ASE read function parameters (ase.io.read), default
       |                   is ':', returning all configurations from the structure.
       |
    Example:
    >>> python ChemicalStructure.py *.xyz ./fake_folder/*.xyz H -ase "{'index' : ':'}" -s
    """


    def __init__(self, chemical_element, xyz_files_paths=None, ase_params={'index' : ':'}):
        self.xyz_files_paths = xyz_files_paths
        self.magres_array_lengths = []
        self.arrays_for_magres = []
        self.structs_list = []
        self.names_list = []
        self.ase_params = ase_params
        self.chemical_element = chemical_element
        self.structs = None
        self.y_pred =None
        self.y_err = None
        self.model = None
        self.csv_file_name = f"{chemical_element}-{(datetime.now().isoformat(sep='_', timespec='seconds'))}-ShiftML2-prediction"
    

    def start_script(self):
        """Start script function 
        
            A function that makes the calls to the appropriate class functions
            to make the chemical shift predictions and output the csv and magers files.

        """

        self.read_xyz_files(self.xyz_files_paths)
        self.load_model()
        self.build_model_representation()
        self.predict_set()
        self.output_csv()
        self.prepare_magres()
        self.output_magres()


    def read_xyz_files(self, list_of_paths):
        """Read xyz file path

            A function that reads all the xyz files in an array of xyz file paths.
            The function creates a list of atom objects collected from all xyz files.
            It also creates an list of file names 'names_list' used to name the magres files.

        | Args:
        |   Input:
        |       list_of_paths (list): list of paths leading to xyz files.
        |
        |   Output:
        |       structs_list (list): list of atom objects found in every xyz file.
        |       names_list (list): List of xyz files of format: 
        |       {xyz_file_name}-{file_index}-atom_obj-{atom_object_index}-ShiftML2-prediction.
        |       magres_array_lengths (list): List of lengths corresponding to
        |                               number of time the chemical element 
        |                               to be predicted, appears in the atom 
        |                               objects found within an xyz file.
        |         
        
        """
        for f_index, xyz_file_path in enumerate(list_of_paths):
            print("Path: ",xyz_file_path)
            try:
                structs = ase.io.read(xyz_file_path, **self.ase_params)
            except:
                raise ValueError(f"Failed to read file. \nWrong 'ase_params' provided: {self.ase_params}")
            if not structs:
                raise ValueError("No structures could be read from the paths provided.")
            file_name = os.path.basename(xyz_file_path)
            index_of_dot = file_name.index('.')
            file_name_without_extension = file_name[:index_of_dot]
            if type(structs) == list:
                for o_index, atom_object in enumerate(structs):
                    self.structs_list.append(atom_object)
                    self.names_list.append(f"{file_name_without_extension}-{f_index:04d}-atom_obj-{o_index:04d}")
                    target_atoms = [atom for atom in atom_object if atom.symbol == self.chemical_element]
                    num_atoms = len(target_atoms)
                    self.magres_array_lengths.append(num_atoms) 
            elif type(structs) == 'ase.atoms.Atoms':
                self.structs_list.append(structs)
                self.names_list.append(f"{file_name}-{f_index:04d}-atom_obj-0000")
                target_atoms = [atom for atom in structs if atom.symbol == self.chemical_element]
                num_atoms = len(target_atoms) 
                self.magres_array_lengths.append(num_atoms)  
            else:
                raise TypeError(f"The type of {structs} should be either 'list' or 'ase.atoms.Atoms'")

 
    def load_model(self):
        """Get model function

            Get ShiftML2 model for specific chemical element.

        | Args:
        |   Input:
        |       chemical_element (str): Element to predict ("H", "C", "N", "O", "S", 
        |                                 "F", "P", "Cl", "Na", "Ca", "Mg" or "K").
        |                    
        |   Output:
        |       model (ShiftML2.ShiftML2.Model): ShiftML2 model class.
        | 

        """
        self.model = sml.Model(self.chemical_element)


    def build_model_representation(self):
        """Build model function

            Build the descriptor representation of atomic environments of a given atom in a set of structures.

        | Args:
        |   Input:
        |       structs_list (list): List of atom objects found in every xyz file path.  
        |                    
        |   Output:
        |       model (ShiftML2.ShiftML2.Model): ShiftML2 model class.
        | 

        """
        self.model.build_representation(self.structs_list)


    def predict_set(self):
        """Real prediction function

            Predict target values on a test set using an ensemble of trained models.

        | Args:
        |   Input:
        |       structs_list (list): List of atom objects found in every xyz file path.  
        |                    
        |   Output:
        |       model (ShiftML2.ShiftML2.Model): ShiftML2 model class.
        |       y_pred2 (numpy.ndarray): Predicted values for the real set.
        |       y_err2 (numpy.ndarray): Predicted errors for the real set.
        | 

        """
        start_time = time.time()
        self.y_pred, self.y_err = self.model.predict()
        end_time = time.time()
        elapsed_time = end_time - start_time
        print("Elapsed time: ", elapsed_time)

            
    def output_csv(self):
        """Output csv function

            Output one csv file for all xyz files found in user provided paths.

        | Args:
        |   Input:
        |       y_pred2 (numpy.ndarray): Predicted values for the real set.
        |       y_err2 (numpy.ndarray): Predicted errors for the real set.
        |
        |   Output:
        |       file (csv): the csv file will be stored within the current directory.
        |                                   
        | 

        """
        prep_array = np.asarray([self.y_pred, self.y_err])
        np.savetxt(f"./{self.csv_file_name}.csv", prep_array, delimiter=",")


    def prepare_magres(self):
        """Prepare magres function

            Slice predicted values returned as one big array into smaller arrays to be
            replaced their corresponding values at corresponding positions within their 
            respective atom object and transform it into a magres file.

        | Args:
        |   Input:
        |      y_pred2 (numpy.ndarray): Predicted values for the real set.  
        |      magres_array_lengths (list): List of lengths corresponding to
        |                               number of time the chemical element 
        |                               to be predicted, appears in the atom 
        |                               objects found within an xyz file.
        |
        |   Output:
        |       array_for_magres (list): List of predicted values for the real set 
        |          corresponding to their atom object found in a specific xyz file.     
        |       arrays_for_magres (list): List of array_for_magres arrays.
        | 

        """
        start = 0
        for slice_length in self.magres_array_lengths:
            end = start + slice_length
            array_for_magres = self.y_pred[start:end]
            self.arrays_for_magres.append(array_for_magres)
            start = end


    def output_magres(self):
        """Output magres function

            Output one magres file for each atom obect found in every xyz file.

        | Args:
        |   Input:
        |       arrays_for_magres (numpy.ndarray): List of lists of predicted values
        |                                          for the real set.  
        |       structs_list (list): List of atom objects found in every xyz file path.
        |       names_list (list): List of xyz files of format: 
        |       {xyz_file_name}-{file_index}-atom_obj-{atom_object_index}-ShiftML2-prediction.
        |       chemical_element (str): Element to predict ("H", "C", "N", "O", "S", 
        |                                 "F", "P", "Cl", "Na", "Ca", "Mg" or "K").
        |       indices (list): label column found in xyz. if label column if missing
        |                 it will be replaced with their corresponding chemical element symbol.
        |
        |   Output:
        |       file (magres): magres files will be stored within the current directory.
        |                                   
        | 

        """
        for i_list, (pred_array, atoms, name) in enumerate(zip(self.arrays_for_magres, self.structs_list, self.names_list)):
            ms = np.zeros((len(atoms),3,3))
            atom_index = 0
            for i_atom, atom in enumerate(atoms):
                if atom.symbol == self.chemical_element:
                    ms[i_atom][np.diag_indices(3)] = pred_array[atom_index]
                    atom_index += 1
            atoms.arrays['ms'] = ms
            if atoms.has('labels'):
                labels = list(atoms.get_array('labels'))
            else:
                labels = atoms.get_chemical_symbols()  
            indices = [labels[:i + 1].count(labels[i]) for i in range(len(labels))]
            atoms.arrays['indices'] = indices
            atoms.info['magres_units'] = {'ms': 'ppm'}
            atoms.write(f'./{name}-ShiftML2-prediction.magres', format='magres')
        print(f"Outputted magres files count: {i_list + 1}")

if __name__ == '__main__':
    cli()
    # sys.exit(cli())


