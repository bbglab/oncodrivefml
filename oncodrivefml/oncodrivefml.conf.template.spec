[signature]
method = option('none', 'full', 'complement', 'bysample', 'file')
classifier = string(default='none')
use_only_mapped_elements = boolean(default=False)
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

[background]
sampling = integer
recurrence = boolean

[statistic]
method = option('amean', 'gmean', 'maxmean', 'amean_scoresmodif')
subs = boolean(default=True)

use_gene_mutations = boolean(default=True)

    [[indels]]
        enabled = boolean(default=False)
        method = option('pattern', 'stop', default='pattern')
        enable_frame = boolean(default=False)

        window_size = integer(default=10)
        weight_function = option('constant', 'linear', 'logistic', default='linear')

        stop_function = option('mean', 'median', 'random', 'random_choice', default='mean')

[settings]
cores = integer(default=None)
drmaa = integer(default=None)
drmaa_maximum = integer(default=100)
queues = list(default=None)
