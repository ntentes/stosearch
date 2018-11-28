---
output:
  html_document: default
  pdf_document: default
---
# stosearch
Stochastic Search Simulation With Resetting

###Description

Allows a user to specify a positive starting point and simulate (data or visually or both) 1 to N searchers moving according to standard brownian motion ($\mu = 0$ and $\sigma^2 = t$) with or without deterministic or stochastic resetting according to rates set by the user. Search will stop when any searcher finds the target (the origin).

####Output

If `output_data==TRUE`, then the function returns a list of named vectors: 

**"TargetFound"** = T or F  
**"TimeUnits"** = time elapsed until a searcher finds the target  
**"SearchCost"** = TimeUnits multiplied by # of searchers  
**"SearchPaths"** = a matrix of numerics where the rows (of length n) are the search path of each searcher  
**"NumResets"** = number of times a searcher was reset (provides NA is no resetting is selected)  
**"ResetTimes"** = a vector of reset times (NA is provided if no resetting is selected). 

If `output_visual==TRUE`, a plot of the search is printed, with each searcher a different colour, a reset marked by a vertical dotted blue line, and a completed search point (the time value where the target is found) is marked by a vertical dashed green line.

###Examples

``` {r}
stochastic_search(x0 = 20, searchers = 5, reset_type = 'deterministic', reset_rate = 1/367, output_data = FALSE)
```

```{r}
stochastic_search(x0 = 20, searchers = 1, t=2000, reset_type = 'stochastic', reset_rate = 1/367, output_data = FALSE)
```

```{r}
stochastic_search(x0 = 100, searchers = 15, reset_type = 'none', output_data = FALSE)
```

```{r}
stochastic_search(x0 = 1, searchers = 3, n=10, reset_type = 'none', output_data = FALSE)
```
