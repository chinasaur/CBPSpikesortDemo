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

(0) Load raw data, stored in an array containing electrode voltages
(each row is a separate electrode).  Sampling rate should be at least
5kHz.

PRE-PROCESSING:

(1) Filter. Purpose is to eliminate low and high frequency noise
(increasing the signal-to-noise ratio of the data) and to allow
isolation of noise regions (for whitening, step 2) by thresholding.

(2) Whiten noise.  The CBP objective function is most easily and
efficiently computed when the background noise is uncorrelated
(white), over both time and electrodes.  Covariance over time and
electrodes is computed on low-amplitude portions of data.  The entire
data array is then (separably) whitened by linearly transforming
across electrodes, and filtering over time.

(3) Select number of cells, and initialize spike waveforms.  Initial
waveforms are obtained from the centroids of k-means clustering in a
principal components space. Note: this NOT used to identify/estimate
spikes - it is only used to determine the number of cells, and to
obtain initial estimates of their waveforms.

CBP SPIKE SORTING:

(4) Partition data into segments (optional).  To improve efficiency
(especially when using multiple cores or machines), the data array is
partitioned into "snippets", separated by spike-free intervals.

(5) Use CBP (see article listed above) to estimate amplitudes/times of
spikes associated with each waveform.

(6) Decide on which spikes to keep/discard, by thresholding the
recovered amplitudes.

(7) Re-estimate waveforms (given the spike times, this is a simple
regression problem).  If changes are substantial, repeat from step (5)
until convergence.

POST-PROCESSING:

(8) Compare CBP results to Clustering, and ground truth (true spikes, if available).

=====================================================================================