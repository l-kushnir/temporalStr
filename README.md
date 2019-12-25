The repository contains Matlab code for simulating four versions of 
efficient balanced spiking network described in the article http://arxiv.org/abs/1912.10262

`FastOnly.m`  - the code for simulating the network with only fast connections and producing Figure 2 of the paper

`SlBalDyn.m` - the code for simulation the network with one type of slow connections and no dimensionality expansion, produces Figure 3 from the paper

`IMDyn.m` - the code for simulating the network with one type of slow connections and 2-fold dimensionality expansion for Figure 4a,b

`DMDyn.m` - the code for simulating the network with two types of slow connections and 3-fold dimensionality expansion, Figure 4c,d

`ActiveF1Fast.m`, `ActiveF1.m`, `ActiveF2.m`, `ActiveF3.m`, `bunch.m`, `dynStim.m` - auxiliary functions

`@EffNetwork`, `@Stimulus` - class folders 
