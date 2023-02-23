import os
import numpy as np
import click

from ase.io import read
from soprano.calculate.nmr import NMRCalculator
from soprano.properties.nmr import MSIsotropy
import matplotlib.pyplot as plt

def validate_input(title, chemical_element, path_to_magres, frequency_range ):
    if not isinstance(title, str):
        raise TypeError(f"The following {title} is not a string.")
    elif len(title) > 30:
        raise ValueError(f"The title provided is too long. Max characters permitted is 30")
    
    if not isinstance(chemical_element, str):
        raise TypeError(f"The following {chemical_element} is not a string. Please provide just one chemical element as string.")
    elif chemical_element not in ["H", "C", "N", "O", "Cl", "Ca", "F", "K", "Mg", "N", "Na", "P", "S"]:
        raise KeyError(f"The following chemical element '{chemical_element}' is not a supported chemical element:[H, C, N, O, Cl, Ca, F, K, Mg, N, Na, P, S].")
    
    if frequency_range:
        try:
            frequency_range = eval(frequency_range)
        except:
            raise TypeError(f"You must provide the dict as a string using double quotes\nTry using: \"{frequency_range}\"") 
        if not isinstance(frequency_range, dict):
            raise TypeError("The frequency range is not a dictionary\nProvide an in double quotes dictionary of the form: \"{'min': 0, 'max': 1000}\"")
    for index, path in enumerate(path_to_magres):
        if not path.endswith(".magres"):
            raise click.BadParameter(f"{path} is not a valid path targetting .magres files.")
        else:
            print(f"Path for plot {index:02d}: {path}")    
    atoms_list = [read(path) for path in path_to_magres]
 
    print(f'Magres files count: {len(path_to_magres)}\n')
    return atoms_list, frequency_range

@click.command()
@click.argument('path_to_magres', nargs=-1, type=click.Path(exists=True))
@click.argument('chemical-element')
@click.option('--freq_broad', '-fb', type=float, required=True, help='The frequency broadening factor to use. [required]')
@click.option('--bins_count', '-b', type=int, default=300, help='The number of bins to use in the spectrum. If not provided, the default is 300.')
@click.option('--reference', '-r', type=str, default=None, help='The reference atom for the chemical shifts.')
@click.option('--title', '-t', type=str, default=None, help='A title for the plot.')
@click.option('--frequency_range', '-fr', default=None, help='A dictionary specifying the minimum and maximum frequency to use in the spectrum. The keys should be "min" and "max".')



def output_spectra(chemical_element,  freq_broad, bins_count, path_to_magres=None, reference=None, title=None, frequency_range=None):
    """
    This function generates a 1D powder NMR spectrum for the given chemical element in the `.magres` files found in the specified path.
    
    Parameters:
    path_to_magres (str, optional): The path to the `.magres` file(s) to use for the calculation.
    chemical_element (str): The chemical element to generate the spectrum for (e.g. "H").
    freq_broad (float): The frequency broadening factor to use.
    bins_count (int): The number of bins to use in the spectrum.
    reference (str, optional): The reference atom for the chemical shifts.
    title (str, optional): A title for the plot.
    frequency_range (dict, optional): A dictionary specifying the minimum and maximum frequency to use in the spectrum. The keys should be "min" and "max".

    Returns:
    None

    Raises:
    ValueError: If the specified path does not lead to a `.magres` file or a folder containing `.magres` files.

    Example:
    >>> python output_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -fr "{'min':0, 'max':1000}"
    """
    atoms_list, frequency_range = validate_input(title = title, chemical_element = chemical_element, path_to_magres = path_to_magres, frequency_range = frequency_range )
    print('Plotting the plots..')
    max_freq = None
    min_freq = None
    for index, (atoms, path) in enumerate(zip(atoms_list, path_to_magres)):
        print(f"Start plotting: plot {index:02d}")
        calc = NMRCalculator(atoms)

        if reference:
            calc.set_reference(reference, chemical_element)
            use_reference = True
        else:
            use_reference = False

        calc.set_powder(N=1)

        if not frequency_range:
            iso = MSIsotropy.get(atoms)
            max_freq=max(iso)*1.1
            min_freq=min(iso)*1.1
            spec_elem_iso, freqs = calc.spectrum_1d(chemical_element, max_freq=max_freq, min_freq=min_freq, bins=bins_count,
                                            freq_broad=freq_broad, use_reference=use_reference)
        else:
            spec_elem_iso, freqs = calc.spectrum_1d(chemical_element, max_freq=frequency_range['max'], min_freq=frequency_range['min'], bins=bins_count,
                                            freq_broad=freq_broad, use_reference=use_reference)
        spec_elem_iso = spec_elem_iso / np.sum(spec_elem_iso)

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(freqs, spec_elem_iso, label='Isotropic shifts')
        if reference:
            ax.set_xlabel(f'{chemical_element} chemical shift (ppm)')
            ax.invert_xaxis()
        elif reference == None:
            ax.set_xlabel(f'{chemical_element} chemical shielding (ppm)')
        ax.set_ylabel('Intensity (a.u.)')
        if title:
            ax.set_title(title)

        if frequency_range and not reference:
            ax.set_xlim(frequency_range['min'], frequency_range['max'])
        elif frequency_range and reference:
            ax.set_xlim(-frequency_range['min'], -frequency_range['max'])
        elif not frequency_range and not reference: 
            ax.set_xlim(min_freq, max_freq)
        elif not frequency_range and reference:
            ax.set_xlim(-min_freq, -max_freq)
        ax.set_xlim()
        ax.set_yticks([])
        ax.legend()

        file_name = os.path.basename(path)
        index_of_dot = file_name.index('.')
        file_name_without_extension = file_name[:index_of_dot]
        plt.savefig(f'./{file_name_without_extension}-spectrum.png', dpi=300)
        plt.close()
        print(' Finished plotting.')
    print("Done!")
if __name__ == '__main__':
    output_spectra()