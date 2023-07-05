# RFlocalfdr

Provides a method for setting the significance level of the MDI (mean decrease in impurity) importances from a random forest model.
Based on an empirical Bayes model. See https://www.biorxiv.org/content/10.1101/2022.04.06.487300v2
Thresholding Gini Variable Importance with a single trained Random Forest: An Empirical Bayes Approach
(Robert Dunne, Roc Reguant, Priya Ramarao-Milne, Piotr Szul, Letitia Sng, Mischa Lundberg, Natalie A. Twine, Denis C. Bauer) for full details.

## Install devtools from CRAN
```r
install.packages("RFlocalfdr")
```

## Or from GitHub:
```r
devtools::install_github("parsifal9/RFlocalfdr", build_vignettes = TRUE)
```

## Usage

```R
library(RFlocalfdr)
vignette("simulated",package="RFlocalfdr")
vignette("Smoking",package="RFlocalfdr")

```

## License

GNU General Public License


