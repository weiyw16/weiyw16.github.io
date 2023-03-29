---
date: 2016-04-09T16:50:16+02:00
title: synthetic 2D elastic walkaway VSP
weight: 20
---
## for Train
<!-- introduction of models -->

| Index of Model | Volume           | Location of Well                         | Location of Shots         | Simulation parameters                       | Notes |
|:---------------|------------------|------------------------------------------|---------------------------|---------------------------------------------|-------|
| [vModel-10]({{%relref "content/seismic_data/synthetic_2d_elastic_norm/fortrain/fortrain_model_10.md" %}})      | 240 x 2048 x 126 | rz=50, dz=2, nz=126; rx=15, dx=50, nx=16 | sz=5, sx=50, ds=50, ns=15 | Ricker wavelet, fm=15 Hz, dt=0.002, nt=3001 | cut   |

## for Test
<!-- introduction of models -->

| Index of Model | Volume           | Location of Well            | Location of Shots         | Simulation parameters                       | Notes |
|:---------------|------------------|-----------------------------|---------------------------|---------------------------------------------|-------|
| vModel-10      | 140 x 2048 x 126 | rz=50, dz=2, nz=126; rx=375 | sz=5, sx=50, ds=5, ns=140 | Ricker wavelet, fm=15 Hz, dt=0.002, nt=3001 | cut   |

