=====================================================================================
||  CBPSpikesortDemo 
=====================================================================================

This package contains matlab code for sorting/estimating spikes of
neurons recorded with one or more extracellular electrodes.  Unlike
commonly used clustering methods, this method can recover temporally
overlapping spikes through use of a sparse inverse algorithm known as
Continuous Basis Pursuit (CBP).  The method is described in this
publication:

   A unified framework and method for automatic neural spike identification.
     C Ekanadham, D Tranchina, and E P Simoncelli. J. Neuroscience Methods,
     vol. 222, pp. 47--55, Jan 2014. DOI: 10.1016/j.jneumeth.2013.10.001
   http://www.cns.nyu.edu/~lcv/pubs/makeAbs.php?loc=Ekanadham13            

For optimization, the code uses the Embedded Conic Solver (ECOS)
package, an interior-point algorithm for second-order cone programming,
written by A. Domahidi, E. Chu, and S. Boyd.  
    See http://web.stanford.edu/~boyd/papers/ecos.html
For convenience, this package is included in the ecos subDirectory.

We also include two example datasets with ground truth (correct spike
times), as obtained from these sites:
*** INCLUDE URLS here ***

TO GET STARTED: we suggest you work through the code in the demonstration script, 
    spikesort_demo/cbp_spikesort_demo_script.m
executing one section at a time.

This directory contains a ChangeLog documenting changes made to the
code, as well as a ToDo list.

* Conceptualization and design by Chaitanya Ekanadham and Eero Simoncelli.
* Original matlab code written by Chaitanya Ekanadham, 12/26/2012.
* Rewritten to use ECOS, and interface updated by Peter H. Li, Fall/Winter 2013.
* Current version available at http://www.cns.nyu.edu/~lcv/software.php

=====================================================================================
OUTLINE OF METHOD:

(0) Load raw data, stored in an array containing electrode voltages over
time.  Sampling rate should be at least 5kHz.

PRE-PROCESSING:

(1) Temporal filtering.  Purpose is to eliminate low and high
frequency noise, increasing the signal-to-noise ratio of the data,
allowing crude initial spike finding by peak thresholding.

(2) Noise whitening.  The CBP objective function is most easily and
efficiently computed when the background noise is uncorrelated
(white), over both time and electrodes.  Covariance over time and
electrodes is computed on low-amplitude portions of data.  The entire
data array is then (separably) whitened by linearly transforming
across electrodes, and filtering over time.

(3) Estimate initial spike waveforms.  For this, we use a simple
k-means clustering method.  Note: this NOT used to identify/estimate
spikes - it is only used for waveform estimation.

(4) Partition data (optional).  To improve efficiency of the method
(esp. when using multiple cores or machines), data array is
partitioned into "snippets", separated by spike-free intervals.

CBP SPIKE SORTING:

(5) Use CBP to estimate spike times associated with each waveform.

(6) Re-estimate waveforms.  If changes are substantial, or if results
show unacceptable refractory violations or notches in
cross-correlograms, repeat from step (5).

POST-PROCESSING:

(7) Compare CBP results to Clustering.  If ground truth (true spikes)
are available, compare to these as well.
=====================================================================================