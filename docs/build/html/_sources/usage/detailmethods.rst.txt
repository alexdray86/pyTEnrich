.. _detailmethods:

Detailed methods
________________


**Binomial Model**

The enrichment is performed using a binomial test. The binomial test is an exact test of the statistical significance of deviations from a theoretically expected distribution, considering two possible outcome. In our case : TE overlap with peak (success) or do not overlap (failure). The statistical hypothesis is :

.. math::
    {\displaystyle H_{0}:p_{obs} =p_{exp}}

where :math:`p_{exp}` is the ratio of genome occupancy of each TE family. :math:`p_{exp}` is therefore the expectation based on genome occupancy, and :math:`p_{obs}` is the observed overlap. Since the probability of having k successes from n trials is given by

.. math::
    P(B = k) = \binom{n}{k} p^k (1 - p)^{n - k}

We can calculate the probability to have at least k success by suming up probabilities, from k success to n success. As often the number of success is on the low edge, we prefer to compute the inverse probability :

.. math::
    P(B >= k) = 1 - P(B < k) = 1 - \sum_{i=0}^{k-1} \binom{n}{i} p^i (1 - p)^{n - i}

This probability is our p-value of having at least k success, given a probability p for the overlap, and n trials. The p-values obtained above are then adjusted with the Benjamin-Hochsberg method to correct for multiple testing.

**Assumptions**

The binomial model itself comes with its associated assumptions :

    - The n trials are mutually independent. The *trials*, in the case of ChIP-seq peaks, as not really independant. Considering the large size of the genome, this assumption seems reasonable.

    - The probability of a given outcome is the same for all n samples. We assume that the probability to touch a given TE family/subfamily is the same for each peak (or input sequence) at each trial. This probability is equal to the ratio of genome occupancy and considered as the probability to touch a TE family/subfamily if we would randomly shuffle sequences. 

    - The only source of variation is simple random and binomial. This assumptions is of course broken as ChIP-seq peaks are by nature not random. It might nevertheless be true for the expected probability if we would random shuffle sequence.

**Two different comparison**

From our perpective, there is a few additional assumptions to be taken into account : 

- We assume that the ratio of genome occupancy corresponds to the expected probability if we would random shuffle sequences. As we are dealing with intervals and not single base pair, the probability might be biased by interval length.

To be as precise as possible, *pyTEnrich* will choose one out of two different comparisons, either *te_in_peak* or *peak_in_te*. 
