---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Commands

All commands can be seen if you run: 
```python
python script_name.py --help 
```
## Visual guide

<iframe width="560" height="315" src="https://www.youtube.com/embed/f1UfwtmJnJY?rel=0&amp;controls=0&amp;showinfo=0" frameborder="0" allowfullscreen></iframe>



## Make predictions
Script name: make_predictions.py

```python
python make_predictions.py *.xyz ./example/ababub_xyz_files/*.xyz H -ase "{'index' : ':'}" -s
```

## Create spectra
Script name: create_spectra.py

```python
python create_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -fr "{'min':10, 'max':35}"
```
## Compare spectra
Script name: compare_spectra.py

```python
python compare_spectra.py *.magres H -t custom_title -fb 0.05 -b 300 -r 0.5 -w [1, 4, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1] -fr "{'min':10, 'max':35}"
```