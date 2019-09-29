# DFTAtom
Density Functional Theory in real space, for atoms

Description will follow soon on https://compphys.go.ro

### PROGRAM IN ACTION

[![Program video](https://img.youtube.com/vi/0wgJyz-M9mI/0.jpg)](https://youtu.be/0wgJyz-M9mI)

### SOME DESCRIPTION

I changed the code to 'shoot' from both 'infinity' and from the nucleus and meet and match the solutions. Still does not work very well for heavy atoms, but for the program purpose it seems to work okish.

As it's LDA and there is the assumption of spherical symmetry (for non occupied shells you don't really have it except for the 'special' cases when you go with LSDA, for the others there is cylindrical symmetry), it kind of works only for noble gases, but you can get quite ok results for other atoms, too.

Here is what I get in a bad scenario, Radon:

```
Step: 30
Energy 1s: -3204.75628989 Num nodes: 0
Energy 2s: -546.577961043 Num nodes: 1
Energy 2p: -527.533025582 Num nodes: 0
Energy 3s: -133.369144877 Num nodes: 2
Energy 3p: -124.172862662 Num nodes: 1
Energy 3d: -106.945006745 Num nodes: 0
Energy 4s: -31.2308037361 Num nodes: 3
Energy 4p: -27.1089853949 Num nodes: 2
Energy 4d: -19.4499946095 Num nodes: 1
Energy 4f: -8.95331838731 Num nodes: 0
Energy 5s: -5.88968288409 Num nodes: 4
Energy 5p: -4.40870276841 Num nodes: 3
Energy 5d: -1.91132962721 Num nodes: 2
Energy 6s: -0.626570717299 Num nodes: 5
Energy 6p: -0.293180028208 Num nodes: 4
Etotal = -21861.3469009 Ekin = 21854.6727109 Ecoul = 8632.01604025 Eenuc = -51966.1203978 Exc = -381.915254228

Finished!

1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6 
```

I used 16 for 'multigrid levels' (that means 65537 nodes) 0.0002 for delta, mixing 0.5 and the max radius 50.
The results are not perfect, I guess with some other parameters they might be improved somewhat. I might improve it more in the future.
Here are the NIST values for comparison: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-radon

For a lighter noble gas I get better results, for Argon for example:

```
Step: 28
Energy 1s: -113.800134596 Num nodes: 0
Energy 2s: -10.794172463 Num nodes: 1
Energy 2p: -8.44343931482 Num nodes: 0
Energy 3s: -0.883384105504 Num nodes: 2
Energy 3p: -0.382330146447 Num nodes: 1
Etotal = -525.946200915 Ekin = 524.969815291 Ecoul = 231.458124025 Eenuc = -1253.13198594 Exc = -29.2421542868

Finished!

1s2 2s2 2p6 3s2 3p6 
```

I used 15 for 'multigrid levels' (that means 32769 nodes), 0.0005 for delta, mixing 0.5 and the max radius 25.
The energy levels results match all given decimals from NIST.

Here is the NIST data for comparison: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-argon

The program is faster than I expected, especially because I went to multi-grid instead of Numerov (or Runge-Kutta) for Poisson, so I expected it to be slower than with the other approach.

