# Bootstrapping pipeline

# Introduction
This is a bootstrapping pipeline.

# Scripts
## wrap_features.R 
The script wraps the input data in a *GenomicRange* object to speed up the calculation of intersections on BED files.

### Inputs:
- *data* (*d*): the given feature file. Currently accepts a Rda file.
- *out* (*o*): The file path of the output file.

# Scripts
## wrap_clusters.R
The script wraps the cluster data in a *GenomicRange* object to speed up the calculation of intersections on BED files.

### Inputs:
- *input* (*i*): The path of the folder containing the cluster data for the chromosomes.
- *out* (*o*): The path of the output folder.
- *pattern* (*p*): The naming patter of the input chromosome files. The pattern should contain a single *%* character which will be parsed to the chromosome name. Example: *%_spec_res.Rda* will expect an input in the shape of *chr1_spec_res.Rda*. 
- *save* (*s*): The naming pattern for the wrapped saved cluster files. Example: a pattern of *%_wrapped.Rda* will save the wrapped data as *chr1_wrapped.Rda*.