
# stosearch

Stochastic Search Simulation With Resetting

![Example of a Stochastic Search](https://i.imgur.com/KfApyEK.png)

### Description

This package can be used for brownian motion search simulations and data can be collected for experimentation. This package can also be used for educational purposes, as it provides good visual cues on the mechanics of a stochastic search with resetting.

This package contains one main function which allows a user to specify a positive starting point and simulate a number of searchers moving according to standard brownian motion (mu = 0 and sigma^2 = t) with or without resetting, until they find the target (which is fixed at the origin).

The function accepts a variety of parameters including the type of resetting (including resetting according to a poisson process or deterministic or none), the number of searchers, the starting point, and others. More parameters will be worked in as the package improves.

Regarding the `reset_rate` parameter, the default is `'optimal'`, which is `1/(2*(x0^2))`, where `x0` is the initial search point. This is proportional to the optimal rate for a stochastic search with resetting as per Bhat, De Bacco & Redner (2016).

The function provides a plot of the search (if wanted) as well as some data collected on the simulation of the search, as explained below.

#### Output

If `output_data==TRUE`, then the function returns a list of named vectors: 

+ **"TargetFound"** = T or F.  
+ **"TimeUnits"** = time elapsed until a searcher finds the target.  
+ **"SearchCost"** = TimeUnits multiplied by # of searchers.  
+ **"SearchPaths"** = a matrix of numerics where the rows (of length n) are the search path of each searcher.  
+ **"NumResets"** = number of times a searcher was reset. For deterministic resetting, each individual reset adds 1 to the count. For deterministic resetting, all searchers are reset at once, so it adds a number to the count equal to the number of searchers (NA is provided if no resetting is selected).  
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


```R
ex5 <- stochastic_search(x0 = 20, searchers = 2, reset_type = 'deterministic', reset_rate = 1/500, output_visual = FALSE)

> ex5$TargetFound
[1] TRUE

> ex5$TimeUnits
[1] 638.77

> ex5$SearchCost
[1] 1277.54

> ex5$NumResets
[1] 2

> ex5$ResetTimes
[1] 500

> dim(ex5$SearchPaths)
[1]      2 100000
```


## Future Improvements

+ Seperate search time vector in search with stochastic resetting so a user can determine which searcher reset at each time  
+ Add the ability to plot an already performed search whose data is saved
+ Other improvements (Feel free to contact me with more!)
