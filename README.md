# Sky
Parametric Modeling Tool for Aircraft Design.

## Description
The program is devoted to BWB design(a.k.a. Blended-Wing-Body), and aims at generating models ready for 
CFD calculation and further optimization.

## Implementation
Generally speaking, the STL model is composed of thousands of triangles, so we just need to craft out
all the geometric elements like curves, lines and faces manually, and the discrete them into triangles.
Once you get all these tiny components, assemble them all and you get the model.

## Progression
2016-10 ~ 2017-01: For the first step, I wrote a simple prototype that generate STL models to get familiar with this project.   
2017-02 ~ Now: Next, I will focus on BWB modeling and using IGES instead of STL.

## Examples
![V1.0模型](https://github.com/cangyu/Sky/blob/master/pic/VIEW_V1.0.JPG)

## Dependencies
If you want to run this program, following packages are needed:
> * VTK
> * PyQt5
> * Numpy+MKL
> * SciPy
> * Matplotlib
> * Numpy-STL
> * Python 3.5 (Suggested)

## Declaration
Just for my postgraduate research.
