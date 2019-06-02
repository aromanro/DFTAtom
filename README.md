# DFTAtom
Density Functional Theory in real space, for atoms

Description will follow soon on https://compphys.go.ro

Warning! It still has issues!

The problem is that the shooting method used is rather crude, just from far away towards the nucleus. I first used a zero approximation at the boundary (justified by the fact that the wavefunction is zero at infinity) and some very small value one step inward. Normalization in the end should take care of adjusting the value to the proper one. A better guess is to use an exponential drop of the wavefunction at large distances from nucleus and that's what is there currently. Unfortunately this is no good for using a large max radius, so for heavy atoms it works up to 10 or so, for big values it does not converge and you get NaNs.

A better guess would be to use WKB approximation and even better would be to use the shooting method from both limits (both inward and outward integration), 'meeting' the solutions at the classical turning point and matching them along with their derivatives.
Another way: Green's function method.

I'm currently quite busy at the moment so what is there should do for now. It's not very hard to improve it, though, so I might do it sometime.
