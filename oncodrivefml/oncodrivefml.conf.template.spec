[genome]
build = option('hg18', 'hg19', 'hg38', default='hg19')

[signature]
method = option('none', 'full', 'complement', 'bysample', 'file', default='full')
classifier = string(default='CANCER_TYPE')
include_mnp = boolean(default=True)
correct_by_sites = option('genome', 'coding', default=None)
use_only_mapped_mutations = boolean(default=False)

path = string(default=None)
column_ref = string(default=None)
column_alt = string(default=None)
column_probability = string(default=None)

[score]
file = string
chr = integer
chr_prefix = string
pos = integer
ref = integer(default=None)
alt = integer(default=None)
score = integer
element = integer(default=None)
extra = integer(default=None)

[statistic]
method = option('amean', 'gmean', default='amean')

sampling = integer(default=100000)
sampling_max = integer(default=1000000)
sampling_chunk = integer(default=100000000)
sampling_min_obs = integer(default=3)

per_sample_analysis = option('amean', 'gmean', 'max', default=None)

subs = boolean(default=True)
mnp = boolean(default=True)

use_gene_mutations = boolean(default=False)

    [[indels]]
        enabled = boolean(default=True)
        method = option(stop', 'max', default='max')

        enable_frame = boolean(default=False)
        max_repeats = integer(default=0)

        stop_function = option('mean', 'median', 'random', 'random_choice', default='mean')

        indels_simulated_with_signature =  boolean(default=False)

[settings]
cores = integer(default=None)
drmaa = integer(default=None)
drmaa_maximum = integer(default=100)
queues = list(default=None)
