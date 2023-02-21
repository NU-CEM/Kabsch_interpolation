# Kabsch interpolation for hybrid materials

<mark>Warning: Always check that the interpolated structures are correct - you can visualise the generated structures using [vesta](https://jp-minerals.org/vesta/en/) or similar.</mark>

If you use this code please consider:
- [Citing the associated paper](https://arxiv.org/abs/2302.08412) (currently under review)
- [Citing the Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/faq.html#how-should-i-cite-ase)
- Letting me know so that I can start a list of where it has been used successfully. I'll link to your project from this page.

If you have any questions or comments about this code please contact l.whalley@northumbria.ac.uk or (preferably) raise an issue on the [repository issue tracker](https://github.com/NU-CEM/Kabsch_interpolation/issues). 

### What does this code do? 🖥️
- This code interpolates between two crystal structures
- It is aimed at hybrid organic-inorganic materials where distortions to the inorganic framework can be accurately described using translational interpolation, whereas molecular species require a description of translation and rotation.
- For isolated atoms it linearly interpolates along the translation vector mapping between the atom start position and atom end position. 
- For molecules it linearly interpolates along the translation vector mapping between the molecule start centre-of-mass and molecule end centre-of-mass, and combines this with a linear interpolation along the rotation vector mapping between the molecule start orientation and molecule end orientation. To do so it uses the [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm). 
- It makes extensive use of the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/index.html) and so can be used to read-in/generate structures for various atomistic simulation codes.

### What does this code assume? ✍🏽
- any molecular distortion can be ignored; the molecules are treated as rigid bodies.
- molecules can be identified using natural cutoffs. Please see the [following link](https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html#ase.neighborlist.natural_cutoffs) for more information. Alternatively, a neighbours list can be explicitly provided.
- atoms in the start Atoms object and end Atoms object are listed in the same order. 

### How do I use this code for my own research? 🔬
- you need to have the [ASE](https://wiki.fysik.dtu.dk/ase/index.html), [scipy](https://scipy.org/) and [numpy](https://numpy.org/doc/stable/index.html) Python libraries installed. If you are new to Python, [Anaconda](https://www.anaconda.com/products/distribution) is recommended.
- you can generate interpolated structures using the [Jupyter Notebook](https://github.com/NU-CEM/Kabsch_interpolation/blob/main/Kabsch_interpolation.ipynb) or the [the script](https://github.com/NU-CEM/Kabsch_interpolation/blob/main/kabsch_interpolation.py). 
- If you are using the script, you can call it from the command line e.g. `python kabsch_interpolation.py POSCAR_start.vasp POSCAR_end.vasp -m "CNH6"`. You can access information about the commnand line arguments using `python kabsch_interpolation.py --help`.

### Something isn't quite right...🙋
- If your molecules are not interpolating as expected you may need to vary the `mic_cutoff` parameter. <mark>Too large or too small values for `mic_cutoff` can lead to unexpected output!</mark>
- ASE can also be used to identify unusual changes in bond lengths. There is example code for this in the final cell of [the notebook](https://github.com/NU-CEM/Kabsch_interpolation/blob/main/Kabsch_interpolation.ipynb).
- If you have any questions about this code please contact l.whalley@northumbria.ac.uk or (preferably) raise an issue on the [repository issue tracker](https://github.com/NU-CEM/Kabsch_interpolation/issues). 
