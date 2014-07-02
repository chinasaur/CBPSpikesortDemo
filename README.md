=====================================================================================
  CBPSpikesortDemo
=====================================================================================

This package contains demonstration code for sorting/estimating spikes
of neurons recorded on one or more extracellular electrodes.  Unlike
commonly used clustering methods, this method can correctly recover
temporally overlapping spikes through use of a sparse inverse
algorithm known as Continuous Basis Pursuit (CBP).  The method is
described in this publication:

   A unified framework and method for automatic neural spike identification.
   C Ekanadham, D Tranchina, and E P Simoncelli. J. Neuroscience Methods,
   vol. 222, pp. 47--55, Jan 2014. DOI: 10.1016/j.jneumeth.2013.10.001
   Available at: http://www.cns.nyu.edu/~lcv/pubs/makeAbs.php?loc=Ekanadham13            

- Conceptualization and algorithm design by Chaitanya Ekanadham and Eero Simoncelli.
- Original matlab code written by Chaitanya Ekanadham, 12/26/2012.
- Rewritten to use ECOS, and interface updated by Peter H. Li, Fall/Winter 2013.
- Current version is available at http://www.cns.nyu.edu/~lcv/software.php

For optimization, the code relies on the Embedded Conic Solver (ECOS)
package, an interior-point algorithm for second-order cone programming,
written by A. Domahidi, E. Chu, and S. Boyd.  See
    http://web.stanford.edu/~boyd/papers/ecos.html
For convenience, this package is included in the ecos subFolder.

We also include two example datasets with ground truth (correct spike
times), as obtained from these sites:
*** INCLUDE URLS ***

=====================================================================================
OUTLINE OF METHOD:

(0) Load raw data, in an array containing electrode voltages over
time.  Sampling rate should be at least 8kHz.

(1) Temporal Filtering.  Purpose is to eliminate low and high
frequency noise, increasing the signal-to-noise ratio of the data, and
allowing crude initial spike finding by thresholding.

(2) Estimate covariance of background noise, and whiten the data
according to this covariance.  The CBP method is based on an objective
function that is most easily and efficiently computed when the
background noise is white.  

(3) Obtain initial estimates of spike waveforms.  For this, we use a
simple k-means clustering method.  Note: this NOT used to
identify/estimate spikes - it is only used for waveform estimation.

(4) Partition data into chunks, separated by spike-free intervals.
This is optional: it is only done to improve efficiency of the
algorithm.

(5) Use CBP to estimate spike times associated with each waveform.

(6) Re-estimate waveforms.  If desired, repeat from step (5). 
=====================================================================================