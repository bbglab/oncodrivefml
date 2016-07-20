The method
==========

OncodriveFML tries to find which gens have mutation with a higher impact
than the expected by random mutations. To do so, the value of the
impact scores of the observed mutations is compared with the value of the
same number of mutations occurring in random positions of the gen.

The method operates differently with substitutions and with insertions or deletions.

Substitutions
-------------

The substitutions that occur in the dataset are grouped by the signature classifier
(see details in :mod:`oncodrivefml.signature`).

Within each group, the different substitutions are counted, taking as identifier the previous base, the
mutated base and the next base of the reference genome and the same triplet but with the altered base.
Those counts are used to compute of probability of each substitution to occur (considering the surrounding bases).


Observed
^^^^^^^^

Each observed substitution has a score obtained from the scores file.

Simulated
^^^^^^^^^

When performing the simulation, the same number of substitutions are simulated.

Any substitution within the same gen can be simulated, but the probability of each comes from the
signature that the observed mutation has.

Finally, the simulation is performed a certain number of times.


Insertions and deletions
------------------------

With insertions and deletions the approach is different,
as they do not have only one change but multiple and their effect can
comprise the rest of the gen.

See :mod:`oncodrivefml.executores.indels` to find details.
