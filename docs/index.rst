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

To compute the statistical enrichment, pyTEnrich compare the overlap between input (typically bed files representing ChIP-seq peaks) and transposable elements subfamilies / families. The goal is the find families with a significant overlap within the input bed files. pyTEnrich uses a binomial model with our two possible outcome being : either TE overlap with peak (success) or no overlap (failure) with a probability p. The number of trials N corresponds to the number of peaks in input bed file. We can then compare the observed overlap to the expected, to compute a fold change and a p-value for the enrichment of each TE subfamilies. 

**Example** 

We consider a small genome with one gene composed of two exons, and two TE families. There is 3 ChIP-seq peaks detected :

.. image:: images/fig1.png

**Step 1 : Compute genome occupancy**

First, pyTEnrich compute genome occupancy for TE families and input bed files. 

.. image:: images/fig2.jpg

Note that TE genome occupancy are pre-computed for the provided TE database. If another TE database is given, or if a genome subset is provided (explained below), it will be re-computed (takes a few minutes). 

