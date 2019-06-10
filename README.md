# DFTAtom
Density Functional Theory in real space, for atoms

Description will follow soon on https://compphys.go.ro

Warning! It still has issues!

### PROGRAM IN ACTION

[![Program video](https://img.youtube.com/vi/0wgJyz-M9mI/0.jpg)](https://youtu.be/0wgJyz-M9mI)

### SOME DESCRIPTION

I changed the code to 'shoot' from both 'infinity' and from the nucleus and meet and match the solutions. Still does not work very well for heavy atoms and for them you have to set a small radius to be able to work, but for the program purpose it seems to work okish.

As it's LDA and there is the assumption of spherical symmetry (for non occupied shells you don't really have it except for the 'special' cases when you go with LSDA, for the others there is cylindrical symmetry), it kind of works only for noble gases, but you can get quite ok results for other atoms, too.

Here is what I get in a bad scenario, Radon:

```
Step: 61
Energy 1s: -3204.75164592 Num nodes: 0
Energy 2s: -546.573211017 Num nodes: 1
Energy 2p: -527.528294579 Num nodes: 0
Energy 3s: -133.364338962 Num nodes: 2
Energy 3p: -124.168060254 Num nodes: 1
Energy 3d: -106.940206512 Num nodes: 0
Energy 4s: -31.2259880433 Num nodes: 3
Energy 4p: -27.1041705521 Num nodes: 2
Energy 4d: -19.4451799426 Num nodes: 1
Energy 4f: -8.94850192618 Num nodes: 0
Energy 5s: -5.88490517484 Num nodes: 4
Energy 5p: -4.4039316565 Num nodes: 3
Energy 5d: -1.9065891349 Num nodes: 2
Energy 6s: -0.622434550978 Num nodes: 5
Energy 6p: -0.28942582178 Num nodes: 4
Etotal = -21861.3437843 Ekin = 21853.9298692 Ecoul = 8632.41882992 Eenuc = -51965.7711944 Exc = -381.921289008

Finished!

1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f14 5s2 5p6 5d10 6s2 6p6
```

I used 12 for 'multigrid levels' (that means 4097 nodes), 0.001 for delta, mixing 0.5 and I have to drop the max radius to 6.
The results are not very good, I guess with some other parameters they might be improved somewhat, but the code is missing something. Here are the NIST values for comparison: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-radon

Later edit: Now the radius can get up to 11 in this case and still converge, with a workaround. I will make this much better when I'll have time, by starting the inwards shooting from the 'proper' position depending on the energy level.

For a lighter noble gas I get better results, for Argon for example:

```
Step: 22
Energy 1s: -113.800135292 Num nodes: 0
Energy 2s: -10.7941726364 Num nodes: 1
Energy 2p: -8.44343975586 Num nodes: 0
Energy 3s: -0.883384108382 Num nodes: 2
Energy 3p: -0.382330185919 Num nodes: 1
Etotal = -525.946165301 Ekin = 524.96616451 Ecoul = 231.458075561 Eenuc = -1253.12827587 Exc = -29.2421295064

Finished!

1s2 2s2 2p6 3s2 3p6 
```

The parameters are as above but max radius is increased to 10. 
Here is the NIST data for comparison: https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations-argon

The program is faster than I expected, especially because I went to multi-grid instead of Numerov (or Runge-Kutta) for Poisson, so I expected it to be slower than with the other approach.


