Output
======

Oncodrivefml outputs a :file:`.tsv` file with the results of the method.

In the file, the following columns can be found:

index
    Gene ID from Ensembl

muts
    number of mutations found in the dataset for that gene

muts_recurrence
    number of mutations that do not occur in the same position

samples_mut
    number of mutated samples in the gene

pvalue
    times that the observed value is higher or equal than the
    expected value, divided by the number of randomizations

qvalue
    ``pvalue`` corrected using the Benjamini/Hochberg correction
    (for samples with at least 2 ``samples_mut``)

pvalue_neg
    times that the observed value is lower or equal than the
    expected value, divided by the number of randomizations

qvalue_neg
    ``pvalue_neg`` corrected using the Benjamini/Hochberg correction
    (for samples with at least 2 ``samples_mut``)

symbol
    HGNC Symbol

