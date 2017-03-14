
Pitfalls
========

Considering that MNP mutations contribute to the
signature as a set of SNPs mutations adds
more mutations to these counts compered
to the ones being analysed because OncodriveFML
simulates a MNP mutation as only one mutation.

If the scores files lacks scores for some positions
or certain alterations, OncodriveFML ignores them.

If, for any reason, your signatures lack certain
triplets (probability 0) that are the only ones present in certain
region, OncodriveFML will not compute a P value
for that region.

OncodriveFML statistical power is limited
by the number of simulations performed in each regions.
You can increase the number of simulations,
but the time cost is exponential.

Indels do not contribute to the signature.
If you choose to simulate indels as substitutions
taken signature into account,
beware that they have not contribute it own signature.

On the other hand, if you choose to no use
signature with the indels, their probability
is the inverse of the number of distinct
trinucleotides for all the regions multiplied
by three (there are 3 possible changes).
Typically the value should be 1/192,
and that is the default value OncodriveFML uses.
If OncodriveFML does signature correction,
it obtains the number of distinct triplets
from the correction.
