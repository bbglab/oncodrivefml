[signature]
method = option('none', 'full', 'complement', 'bysample', 'file')
classifier = string(default='none')
path = string(default=None)
column_ref = string(default=None)
column_alt = string(default=None)
column_probability = string(default=None)

[score]
file = string
chr = integer
chr_prefix = string
pos = integer
ref = integer
alt = integer
score = integer
element = integer(default=None)
extra = integer(default=None)

[background]
sampling = integer
recurrence = boolean
range = integer(default=None)

[statistic]
method = option('amean', 'gmean', 'maxmean')
subs = option('enabled', 'none', default = 'enabled')
indels = option('max', 'none', default = 'none')
indels_max_repeats = integer(default=3)

[settings]
cores = integer(default=None)
drmaa = integer(default=None)
drmaa_maximum = integer(default=100)
queues = list(default=None)
