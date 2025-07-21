This repository containts functions for estimating and conducting (bootstrapped) inference of multiple causal effects from a single source of exogenous variation. The method utilizes a 'two stage residual inclusion' (2SRI) style estimator paired with flexible first stage fits, exploiting non-linearities in the relationship between the source of exogenous variation and the endogenous variables to enhance identification. A diagnostic function is provided for evaluating whether the necessasry rank condition has been satisfied.

To import this code, use

```R
source("https://raw.githubusercontent.com/christophercschwarz/sharpeningbluntinstruments/refs/heads/main/sharpening_blunt_instruments_function.R")
```
