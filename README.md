# Paper_code
## code for paper

### Requried R packages
```
data.table >= 1.12.8
intervals >= 0.15.2
```

### Code explanation 
##### explanation the function of ker_ts_fixed.R code
```
ker_ts_fixed.R produces Figure 5 in the Web Appendices, and Figure 4 and 6 in Respose to the comments
```
```
Illustrate the main functions contained in ker_ts_fixed.R code, as follows:

The defined "data.gen()" function is to generate the simulation data set where the obaservation time is fixed;

The defined "beta.fun.ker()" and "beta.fun.ts()" functions are the Newton-Raphson iteration algorithm for the proposed method and the two-stage approach, respectively;

The defined "shape.fun()" function is the estimation of the shape function which is useful to estiamte the baseline function in the two-stage approach;

The defined "cvscore.ker()" and "cvscore.ts()" functions are the cross-validation score function of the proposed method and the two-stage approach, respectively;

The defined "choose.hb.ker()" and "choose.hb.ts()" functions are the cross-validation approach for bandwidth selection of the proposed method and the two-stage approach, respectively;

The defined "simulation.ker()" and "simulation.ts()" are the simulation for 500 replications of the proposed method and the two-satge approach, respectively;

The ploted figure saved as "cpfixed.eps" produces Figure 5 in the Web Appendices, and Figure 4in Respose to the comments;

The ploted figure saved as "fixkxq.eps" produces Figure 6 in Respose to the comments.
```
##### explanation the function of ker_ts_fixed.R code
```
The defined function "data.gen" 
```