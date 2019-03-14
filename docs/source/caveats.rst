
Caveats
=======

Signature computation is performed using all mutations
in your input file, not only the ones
that map to the region of interest.

If the scores files lacks scores for some positions
or certain alterations, OncodriveFML ignores them.

If, for any reason, your signatures lack certain
triplets (probability equal to 0) that are the only ones present in certain
region, OncodriveFML will not compute a P-value
for that region.

OncodriveFML statistical power is limited
by the number of simulations performed in each regions.
You can increase the number of simulations,
but be aware that the time cost is exponential.

Indels do not contribute to the signatures.
You can simulate indels as substitutions and perform the 
simulations taking the signatures into account, but
be aware that the signatures are not calculated considering indels.

Depending on the values of ``sampling_min_obs`` and
``sampling_chunk``  in the configuration file
the number of simulations performed
for a particular genomic element can differ.
