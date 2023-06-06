
cdef extern from "stdlib.h":
    double drand48() nogil
    void srand48(long int seedval) nogil


def walker_sampling(long samples, long muts, double obs_val, double [:] scores, double [:] probs, long [:] inx, long seed):
    cdef long obs=0, neg_obs=0, i=0, j, size=len(scores)
    cdef double mean, u
    cdef list mean_values = []

    srand48(seed)

    while i < samples:
        mean = 0.0

        for r in range(muts):
            u = drand48()
            j = <long> (drand48() * size)

            if u <= probs[j]:
                mean += scores[j]
            else:
                mean += scores[inx[j]]

        mean = mean / muts
        mean_values.append(mean)

        if mean >= obs_val:
            obs += 1
        if mean <= obs_val:
            neg_obs += 1

        i = i + 1
    
    return obs, neg_obs, mean_values