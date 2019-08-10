# BP-reconstruction
When the number of observations is insufficient, we use compressed sensing reconstruction method to reconstruct high-dimensional signals.
All codes can run in Matlab.

When there is no noise, the data are restored with l1eq_pd.m and BP_no_noise.m files.

// Operating method

Signal length, number of spikes in the signal,  number of observations can be configured. Run BP_no_noise.m directly. Or you could run BPnonoise_jiaoben.m to gain the multiple runs results, and data stored in Excel sheet.

There are two ways to run l1eq_pd.m. The default is the first normal mode, and the second is large scale. It should be noted that when using largescale, the cgsolve.m function in the toolbox needs to be used extra. 
////////////////////////////////////////////////////////////////////////

Three files, BPDN, m, l1qc_logbarrier.m and l1qc_newton.m, are used to recover data in case of noise.

// Operating method

Signal length, number of spikes in the signal,  number of observations can be configured. Run BPDN.m directly. Or you could run BPNDjiaoben to gain the multiple runs results, and data stored in Excel sheet.



The l1qc_logbarrier. m has two invocation modes. The default is the first normal mode and the second is large scale. It should be noted that when using largescale, the cgsolve. m function in the toolbox is also needed.
