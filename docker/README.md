# Using Docker with this project

Two Docker images are available from `envest/rnaseq_titration_results`.
They are highly similar but are based on different versions of R.

### R-4.1.2 version

We recommend using this version for running any analysis.
This version is maintained and will be the one to get updated in the future.

To pull this image, use the tag `R-4.1.2`:

```
docker pull envest/rnaseq_titration_results:R-4.1.2
```

### R-3.6.3 version

We also have an image based on R version 3.6.3.
This image is more representative of the development environment used in earlier (pre-2022) iterations of this analysis and we retain it for posterity.

:warning: We do _not_ recommend using this version for running analysis since recent code updates have changed some behaviors under older versions of R and corresponding package versions. 

To pull this image, use the tag `R-3.6.3`:

```
docker pull envest/rnaseq_titration_results:R-3.6.3
```
