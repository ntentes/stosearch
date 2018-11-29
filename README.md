
# stosearch

Stochastic Search Simulation With Resetting

### Description

This package can be used for brownian motion search simulations and data can be collected for experimentation. This package can also be used for educational purposes, as it provides good visual cues on the mechanics of a stochastic search with resetting.

This package contains one main function which allows a user to specify a positive starting point and simulate a number of searchers moving according to standard brownian motion (mu = 0 and sigma^2 = t) with or without resetting, until they find the target (which is fixed at the origin).

The function accepts a variety of parameters including the type of resetting (including resetting according to a poisson process or deterministic or none), the number of searchers, the starting point, and others. More parameters will be worked in as the package improves.

The function provides a plot of the search (if wanted) as well as some data collected on the simulation of the search, as explained below.

#### Output

If `output_data==TRUE`, then the function returns a list of named vectors: 

+ **"TargetFound"** = T or F  
+ **"TimeUnits"** = time elapsed until a searcher finds the target  
+ **"SearchCost"** = TimeUnits multiplied by # of searchers  
+ **"SearchPaths"** = a matrix of numerics where the rows (of length n) are the search path of each searcher  
+ **"NumResets"** = number of times a searcher was reset (provides NA is no resetting is selected)  
+ **"ResetTimes"** = a vector of reset times (NA is provided if no resetting is selected). 

If `output_visual==TRUE`, a plot of the search is printed, with each searcher a different colour, a reset marked by a vertical dotted blue line, and a completed search point (the time value where the target is found) is marked by a vertical dashed green line. The target (the origin) is always plotted and is marked by the dashed red line.

### Examples

```R
stochastic_search(x0 = 20, searchers = 5, reset_type = 'deterministic', reset_rate = 1/367, output_data = FALSE)
```

![Example 1](https://i.imgur.com/a7X5yNR.png)

```R
stochastic_search(x0 = 20, searchers = 1, t=2000, reset_type = 'stochastic', reset_rate = 1/367, output_data = FALSE)
```

![Example 2](https://i.imgur.com/cM47dB9.png)

```R
stochastic_search(x0 = 100, searchers = 15, reset_type = 'none', output_data = FALSE)
```

![Example 3](https://i.imgur.com/dYEh5vD.png)

```R
stochastic_search(x0 = 20, searchers = 3, n=10, reset_type = 'none', output_visual = FALSE)
```

![Example 4](https://i.imgur.com/OrRIrd1.png)
