# OpenSEES-BoucWenGG
Implementation in OpenSEES of the Bouc-Wen model modified by Gerolymos and Gazetas

Reference: Gerolymos N. and Gazetas G. **Constitutive model for 1-D cyclic soil behaviour applied to seismic analysis of layered deposits**. _Soils and Foudations_, 45-147 (2005).


## Usage (TLC):
Include [`BoucWenGG.dll`](/releases/latest) file inside the same folder of the TCL model file to use the material. The material can be defined by the TCL command:
```tcl
uniaxialMaterial BoucWenGG $matTag $alpha $k0 $strainY $n $gamma $beta $s1 $s2 $mkur
```
| parameter | description |
| --- | --- |
| `matTag` | integer tag identifying the material |
| `alpha` | ratio of the post-yielding stiffness to the initial elastic stiffness |
| `k0` | initial elastic stiffness |
| `strainY` | yielding strain $\varepsilon_y$ |
| `n` | parameter that controls the transition between the initial and post-yielding branch. As $n$ increases the transition becomes sharper (when $n\to\infty$ the model is piece-wise linear). |
| `gamma`, `beta` | parameters that control the shape of the hysteresis loop |
| `s1`, `s2` | control the stiffness degradation upon reversal |
| `mkur` | account for stress and tangent stiffness modification proposed by Drosos et al. (2012). `mkur`==1 modification is aplied, `mkur`==0 modification is not applied |
