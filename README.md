# OpenSEES-BoucWenGG
Implementation in OpenSEES of the Bouc-Wen model modified by Gerolymos and Gazetas (2005).

Author: _Andrea Marchi_ (andrea.marchi@uniroma1.it)

References:
- Gerolymos N. and Gazetas G. **Constitutive model for 1-D cyclic soil behaviour applied to seismic analysis of layered deposits**. _Soils and Foudations_, 45-147 (2005).
- Marchi A. **Improved Bouc-Wen model implementation in OpenSees**. _2nd Eurasian Conference on OpenSees (OpenSees Days)_, Turin (2022).


## Usage (TLC):
Include [`BoucWenGG.dll`](https://github.com/mrc-tech/OpenSEES-BoucWenGG/releases/latest/download/BoucWenGG.dll) file inside the same folder of the TCL model file to use the material. The material can be defined by the TCL command:
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
| `mkur` | account for stress and tangent stiffness modification proposed by Drosos et al. (2012). `mkur`=1 modification is aplied, `mkur`=0 modification is not applied |

## Main behaviour
$F(u(t)) = \alpha \ k_0 \ u(t) + (1-\alpha) \ k_0 \ u_y \ \zeta(t)$

$\dot{\zeta}(t) = \frac{\dot{u}(t)}{u_y} [(\gamma+\beta) - |\zeta(t)|^n \ (\gamma + \beta \ sgn(\dot{u}(t) \ \zeta(t))) ] $


![immagine](https://github.com/mrc-tech/OpenSEES-BoucWenGG/assets/74192712/d8e4b931-d50b-4702-9d57-95e391955b59)
