# rnaseq
this is my Q&D rna-seq pipeline using salmon at its core.

it uses salmon to provide fast quantitation of a couple flavors of reference
indices to generate meaningful QC metrics for rna-seq experiments. it passes
read data through tee'ed FIFO files so it is fully compatible with streaming
input (e.g. downloading rna-seq directly from a data repository, so you don't
need to store it locally).
