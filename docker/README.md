# Dockerfiles

When updating code to run in 2021, we experienced different (worse) prediction behavior than expected based on previous results.

### Issue
See: https://github.com/greenelab/RNAseq_titration_results/issues/41

### Game plan
For now, we need to use a Dockerfile that builds an image with R-3.6.3 (`docker/R-3.6.3/Dockerfile_R-3.6.3`).
In the future, we would like to understand why builds with R-4 failed and use the most up-to-date version of R.
The R-3.6.3 Dockerfile will be developed independently from the R-4 version.
The R-4 version will become outdated and need to be brought back up to speed when the time comes.
