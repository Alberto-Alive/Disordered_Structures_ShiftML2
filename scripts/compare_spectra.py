import numpy as np
import seaborn as sns

from ase.io import read
from soprano.properties.nmr import MSIsotropy
from soprano.calculate.nmr import NMRCalculator

import matplotlib.pyplot as plt
import click

def validate_input(title, chemical_element, path_to_magres, frequency_range, weights):
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

    if weights:
        try:
            weights = eval(weights)
        except:
            raise TypeError(f"You must provide the list as a string using double quotes\nTry using: \"{weights}\"")
        if not isinstance(weights, list):
            raise TypeError(f"Provided {weights} is not a list, use format: -w [0,1,1]")
        elif weights and (len(path_to_magres) != len(weights) or min(weights) < 0):
            raise ValueError(f"You need an arrray of weights (positive values including 0) that is of length {len(path_to_magres)} \nYour array of weights has a length of {len(weights)}:{weights}") 
    for index, path in enumerate(path_to_magres):
        if not path.endswith(".magres"):
            raise click.BadParameter(f"{path} is not a valid path targetting .magres files.")
        else:
            print(f"Path for plot {index:02d}: {path}")    
    atoms_list = [read(path) for path in path_to_magres]
    
    print(f'Magres files count: {len(path_to_magres)}\n')
    return atoms_list, frequency_range, weights
   
@click.command()
@click.argument('path_to_magres', nargs=-1, type=click.Path(exists=True))
@click.argument('chemical-element')
@click.option('--title', '-t', type=str, required=True, help='Title of the average spectra. [required]')
@click.option('--freq-broad', '-fb', type=float, default=0.05, help='Broadening applied to the spectra. If not provided, the default is 0.05.')
@click.option('--bins-count', '-b', type=int, default=300, help='Number of bins in the spectra. If not provided, the default is 300.')
@click.option('--weights', '-w', default=None, help='List of weights (N=[0, ∞]) for each spectra. If not provided, a mean is taken. Weights are not normalised, for example given the weights [0, 1, 1] spectra of index 0 will not count but spectra of index 1 and 2 will count 50% each.')
@click.option('--reference', '-r',  type=float, default=None, help='Reference to use for shifting the spectra. If not provided, no reference is used.')
@click.option('--frequency-range', '-fr', default=None, help='Dictionary specifying the maximum and minimum frequency range to use (use double quotes around the dict). If not provided, range is calculated using chemical element.')


def compare_spectra(chemical_element, title, freq_broad, bins_count, weights=None, path_to_magres=None, reference=None, frequency_range=None):
    """
    This function computes the average spectra of chemical element over multiple ".magres" files.
    
    Parameters:
    path_to_magres (str, required): Path to the folder containing the ".magres" files. Defaults to None.
    chemical_element (str, required): Chemical element to analyze.
    title (str): Title of the average spectra.
    freq_broad (float): Broadening applied to the spectra.
    bins_count (int): Number of bins in the spectra.
    weights (list, optional): List of weights for each spectra. If not provided, a mean is taken. Defaults to None. Weights are not normalised, for example given the weights [0, 1, 1] spectra of index 0 will not count but spectra of index 1 and 2 will count 50% each.
    reference (float, optional): Reference to use for shifting the spectra. If not provided, no reference is used. Defaults to None.
    frequency_range (dict, optional): Dictionary specifying the maximum and minimum frequency range to use. If not provided, range is calculated using chemical element. Defaults to None.
    
    Returns:
    None: The function saves the subplots of individual spectra and the average spectra as PNG files in the current directory.
    
    Raises:
    ValueError: If the path to the ".magres" files is not provided, or if the path leads to less than 2 or more than 30 ".magres" files.
    ValueError: If the length of weights array is not equal to the number of spectra, or if there are negative values in the weights array.
    
    Example:
    >>> python compare_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -w [0,1,1] -fr "{'min':0, 'max':1000}"
    >>> compare_spectra("title", "C", 0.05, 300, reference="N", frequency_range={"min":0, "max":1000})
    """
    atoms_list, frequency_range, weights = validate_input(title = title, chemical_element = chemical_element, path_to_magres = path_to_magres, frequency_range = frequency_range, weights = weights)
    
    #Plotting subplots
    print('Plotting the subplots..')
    iso_spectras = []
    freq_spectras = []
    max_freq = None
    min_freq = None
    for atoms in atoms_list:
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
        iso_spectras.append(spec_elem_iso)
        freq_spectras.append(freqs)

    number_of_subplots = len(iso_spectras)
    fig, axs = plt.subplots(nrows=number_of_subplots, sharex=True, figsize=(12, 22))
    plt.subplots_adjust(hspace=1.5)
    fig.suptitle("Subplots of magres", fontsize=18, y=0.95)

    for elem_freq, elem_iso, ax in zip(freq_spectras, iso_spectras, axs.ravel()):
        ax.plot(elem_freq, elem_iso)
        ax.set_ylabel('Intensity (a.u.)')
    
    if reference:
        axs[-1].set_xlabel(f'{chemical_element} chemical shift (ppm)')
        ax.invert_xaxis()
    elif reference == None:
        axs[-1].set_xlabel(f'{chemical_element} chemical shield (ppm)')
    if frequency_range and not reference:
        ax.set_xlim(frequency_range['min'], frequency_range['max'])
    elif frequency_range and reference:
        ax.set_xlim(-frequency_range['min'], -frequency_range['max'])
    elif not frequency_range and not reference: 
        ax.set_xlim(min_freq, max_freq)
    elif not frequency_range and reference:
        ax.set_xlim(-min_freq, -max_freq)
    plt.savefig(f'./{title}-{number_of_subplots}-spectrum_subplots.png', dpi=300)
    plt.close()
    print(' Finished plotting.')

    #Averaging plot
    print('Averaging subplots..')
    if weights:
        average_iso = np.average(iso_spectras, weights=weights, axis=0)
    else:
        average_iso = np.mean(iso_spectras, axis=0)
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(freqs, average_iso)

    if frequency_range and not reference:
        ax.set_xlim(frequency_range['min'], frequency_range['max'])
    elif frequency_range and reference:
        ax.set_xlim(-frequency_range['min'], -frequency_range['max'])
    elif not frequency_range and not reference: 
        ax.set_xlim(min_freq, max_freq)
    elif not frequency_range and reference:
        ax.set_xlim(-min_freq, -max_freq)

    if reference:
        ax.set_xlabel(f'{chemical_element} chemical shift (ppm)')
        ax.invert_xaxis()
    elif reference == None:
        ax.set_xlabel(f'{chemical_element} chemical shielding (ppm)')
    
    ax.set_ylabel('Intensity (a.u.)')
    fig.suptitle("Average of magres subplots", fontsize=18, y=0.95)
    plt.savefig(f'./{title}-average-spectrum.png', dpi=300)
    plt.close()
    print(' Finished averaging.')

    #Correlation Matrix
    print('Calculating correlation..')
    corr_matrix = np.corrcoef(iso_spectras)
    n = corr_matrix.shape[0]
    groups = []
    for i in range(n):
        for j in range(i+1, n):
            if corr_matrix[i, j] >= 0.99:
                groups.append((i, j))

    unmatched = list(range(n))
    all_groups = []
    while len(unmatched) > 0:
        i = unmatched.pop(0)
        matches = [i]
        inds = [j for j in unmatched if corr_matrix[i][j] >= 0.99]
        matches += inds
        unmatched = [u for u in unmatched if u not in matches]
        all_groups.append(matches)

    plot_title = ''
    for group in all_groups:
        plot_title += ' ' + str(group)

    plt.figure(figsize=(8, 6))
    sns.heatmap(corr_matrix, annot=True, cmap="YlGnBu", annot_kws={"size": 8})
    plt.xlabel("Spectrum index")
    plt.ylabel("Spectrum index")
    plt.title(plot_title)
    plt.savefig(f'./{title}-correlation-matrix.png', dpi=300)
    plt.close()
    for i, fname in enumerate(path_to_magres):
        print(i, fname)
    print(' Finished calculation.')

if __name__ == '__main__':
    compare_spectra()