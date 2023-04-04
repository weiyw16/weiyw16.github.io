---
title : "Pre-processing"
date : 2023-04-04T18:28:49+08:00
---
## Example
Four basic pre-processing tools are used when constructing norm synthetic seismic data from raw simulated data. These four steps are:

* re-sample
* time power
* window
* normalization

In a partial Madagascar `Sconstruct` example as what shown in the following codeblock, `sfbandpass`, `sfwindow` are used for re-sampling, where `sfbandpass` are applied before down-sampling to avoid aliasing. 

Then, `sfpow` is used for amplitude balance along the `time` axis, where attenuation caused by geometry spreading would be alleviated. 

Afterwards, a time window is used to cut the data into the same size, which named as the `height` in common AI applications. Here we have the height 2048 points. 

Finally, we normalize data by traces via `sfscale` to limit the values of all seismic data within negative to positive one.

```python
# parameters
dt = 0.0005 # unit: s
fa = 0.1 # time point of first arrival on one trace 
nt = 2048 # length of time axis
path = os.path.join(cwd, pown_dir) # local path to save files

# Madagascar function
Flow(out_name, in_name, '''
                 put o1=500 d1=20 label1=Depth unit1=m
                 o2=0 d2=%g label2=Time unit2=s |
                 transp plane=12 |
                 bandpass fhi=500 |
                 window j1=4 |
                 pow pow1=2 |
                 window f1=%g n1=%g |
                 scale axis=1
                 datapath=%s
                 ''' % (dt, int(fa/dt/4), nt, path) 
    )
```

A complete workflow of constructing norm data from raw synthetic seismic data can be found at [scripts](https://github.com/weiyw16/weiyw16.github.io/tree/main/scripts/workflow/2-template_produce_data).