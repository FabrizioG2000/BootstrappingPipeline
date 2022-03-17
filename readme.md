# Bootstrapping pipeline

# Introduction
This is a bootstrapping pipeline.

# Scripts
## wrap_features.R 
The script wraps the input data in a *GenomicRange* object to speed up the calculation of intersections on BED files.

### Inputs:
- *data* (*d*): the given feature file. Currently accepts a Rda file.
- *out* (*o*): The file path of the output file.