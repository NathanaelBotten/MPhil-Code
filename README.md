# MPhil-Code
Code used to complete the research I've done for my masters of Philosophy at Adelaide University

Included in this repository is:

1. The code I used to fit the coupling constants of the QMC model to the binding of Lambda hyperons, as in chapter 3 of my thesis. This was done using the Numerov algorithm, and done in Jupyter notebook.
2. The code I used to calculate the binding energy of Cascade hyperons, as in chapter 4 of my thesis. This was also done using the Numerov algorithm, and in a Jupyter notebook.
3. Code to calculate the mean value of the sigma meson field, when considering a Sigma hyperon bound in Helium nucleus. Once again using the Numerov algorithm and done using a Jupyter notebook.
4. Fortran code to calculate the shift in the energy levels in Iron and Carbon nuclei, when considering the binding of a cascade nucleus.

Not included is the values calculated for the finite Coulomb potential in the iron and carbon nuclei. There are however instructions for how to calculate these included in the findeigenvalue.f90 files for both the carbon and iron nuclei. It is recommended to write these to a file, as these are rather computationally intensive processes, that are relatively lengthy.
