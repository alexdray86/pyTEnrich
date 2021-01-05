.. pyTEnrich documentation master file, created by
   sphinx-quickstart on Wed Dec 30 11:29:02 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyTEnrich's documentation!
=====================================

A code for to compute statistical enrichment of transposable elements
on a group of bed files. Code written by Alexandre Coudray from the
laboratory of Virology and Genetics at the EPFL in 2019.

.. toctree::
   :maxdepth: 2
   :caption: Contents: 
   
   usage/installation.rst
   usage/execution.rst
   usage/otherspecies.rst

Overview of the methods
=======================

To compute the statistical enrichment, pyTEnrich compare the overlap between input bed files and transposable elements families. pyTEnrich uses a binomial model to compare the observed overlap to the expected one, given a probability computed from the genome occupancy of TE families. The goal is to compute a fold change and a p-value for the overlap observed between each TE families and a group of input bed files. 

**Example** 

We consider a small genome with one gene composed of two exons, and two TE families. There is 3 ChIP-seq peaks detected :

.. image:: images/fig1.png

**Step 1 : Compute genome occupancy**

First, pyTEnrich compute genome occupancy for TE families and input bed files. 

.. image:: images/fig2.jpg

Note that TE genome occupancy are pre-computed for the provided TE database. If another TE database is given, or if a genome subset is provided (explained below), it will be re-computed (takes a few minutes). 

**Step 2 : Intersect TE and Input bed files and count overlap**

Using Bedtools intersect (link to website), we compute a stringent overlap between input bed files and TE database. The observed overlap are counted to be compared with the expectations.

.. image:: images/fig3.jpg

**Step 3 : Compute the enrichment of TE subfamily / family**

For each input bed files, and enrichment is performed using a binomial test. The two possible outcome according to our model is : either TE overlap with peak (success) or no overlap (failure) with a probability p. p was calculated using genome occupancy of each TE subfam. The number of trials N corresponds to the number of peaks in input bed file. 

As the binomial probability of having k successes from n trials is given by :

.. math::
    P(B = k) = \binom{n}{k} p^k (1 - p)^{n - k}

We can calculate the probability to have at least k successes by suming up probabilities :

.. math::
    P(B >= k) = \sum_{i=k}^n \binom{n}{i} p^i (1 - p)^{n - i}

As often the number of success is on the low edge, we prefer to compute the inverse probability :

.. math::
    P(B < k) = 1 - P(B >= k) = 1 - \sum_{i=0}^{k-1} \binom{n}{i} p^i (1 - p)^{n - i}

