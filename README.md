# DFTAtom
Density Functional Theory in real space, for atoms

Description will follow soon on https://compphys.go.ro

Warning! It still has issues!

### PROGRAM IN ACTION

[![Program video](https://img.youtube.com/vi/0wgJyz-M9mI/0.jpg)](https://youtu.be/0wgJyz-M9mI)

### SOME DESCRIPTION

I changed the code to 'shoot' from both 'infinity' and from the nucleus and meet and match the solutions. Still does not work very well for heavy atoms, but for the program purpose it seems to work okish.

As it's LDA and there is the assumption of spherical symmetry (for non occupied shells you don't really have it except for the 'special' cases when you go with LSDA, for the others there is cylindrical symmetry), it kind of works only for noble gases, but you can get quite ok results for other atoms, too.

Here is what I get in a bad scenario, Radon:

```
Step: 51
Energy 1s: -3204.75868481 Num nodes: 0
Energy 2s: -546.578874176 Num nodes: 1
Energy 2p: -527.534187689 Num nodes: 0
Energy 3s: -133.369443983 Num nodes: 2
Energy 3p: -124.173207469 Num nodes: 1
Energy 3d: -106.945348538 Num nodes: 0
Energy 4s: -31.2309033805 Num nodes: 3
Energy 4p: -27.109092456 Num nodes: 2
Energy 4d: -19.4500920358 Num nodes: 1
Energy 4f: -8.95339416469 Num nodes: 0
Energy 5s: -5.88971691991 Num nodes: 4
Energy 5p: -4.40873567403 Num nodes: 3
Energy 5d: -1.91135386428 Num nodes: 2
Energy 6s: -0.626576190101 Num nodes: 5
Energy 6p: -0.293182968877 Num nodes: 4
Etotal = -21861.3362922 Ekin = 21849.8454244 Ecoul = 8631.97981041 Eenuc = -51961.2538851 Exc = -381.907641938

Finished!

1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6 
```

I used 12 for 'multigrid levels' (that means 4097 nodes), 0.001 for delta, mixing 0.5 and the max radius 15.
The results are not very good, I guess with some other parameters they might be improved somewhat, but the code is missing something. Here are the NIST values for comparison: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-radon

For a lighter noble gas I get better results, for Argon for example:

```
Step: 22
Energy 1s: -113.800145224 Num nodes: 0
Energy 2s: -10.7941743912 Num nodes: 1
Energy 2p: -8.44344188769 Num nodes: 0
Energy 3s: -0.883384239216 Num nodes: 2
Energy 3p: -0.382330261694 Num nodes: 1
Etotal = -525.946124261 Ekin = 524.961376334 Ecoul = 231.457986311 Eenuc = -1253.12339158 Exc = -29.2420953215

Finished!

1s2 2s2 2p6 3s2 3p6 
```

The parameters are as above. 
Here is the NIST data for comparison: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-argon

The program is faster than I expected, especially because I went to multi-grid instead of Numerov (or Runge-Kutta) for Poisson, so I expected it to be slower than with the other approach.


