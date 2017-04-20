
Caveats
=======

MNP mutations contribute to the
signatures as a set of independent SNPs mutations.
This means that the calculation of the signatures is
made with a higher number of mutations compared
to the observed substitutions being analysed because OncodriveFML
simulates a MNP mutation as a single SNP mutation.

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

On the other hand, if you choose to not use the
signatures with the indels, their probability
is the inverse of the number of distinct
trinucleotides for all the regions multiplied
by three (there are 3 possible changes).
Typically the value should be 1/192,
and that is the default value OncodriveFML uses.
If OncodriveFML corrects the signatures,
it obtains the number of distinct triplets
from the correction.

Depending on the values of ``sampling_min_obs`` and
``sampling_chunk``  in the configuration file
the number of simulations performed
for a particular genomic element can differ.
