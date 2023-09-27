# OpenSEES-BoucWenGG
Implementation in OpenSEES of the Bouc-Wen model modified by Gerolymos and Gazetas

Reference: Gerolymos N. and Gazetas G. **Constitutive model for 1-D cyclic soil behaviour applied to seismic analysis of layered deposits**. _Soils and Foudations_, 45-147 (2005).


## Usage (TLC):
```tcl
uniaxialMaterial BoucWenGG $matTag $alpha $k0 $strainY $n $gamma $beta $s1 $s2 $mkur
```
| `matTag` | integer tag identifying the material |
| `alpha` | the |
