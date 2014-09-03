#######################################################################
#               _______     ______       _     _________              #
#               |_   __ \ .' ___  |     / \   |  _   _  |             #
#               | |__) |/ .'   \_|     / _ \  |_/ | | \_|             #
#               |  ___/ | |           / ___ \     | |                 #
#               _| |_    \ `.___.'\ _/ /   \ \_  _| |_                #
#              |_____|    `.____ .'|____| |____||_____|               #
#                                                                     #
#                                                                     #
#                       _,'|             _.-''``-...___..--';)        #
#                      /_ \'.      __..-' ,      ,--...--'''          #
#                    <\    .`--'''       `     /'                     #
#                      `-';'               ;   ; ;                    #
#                __...--''     ___...--_..'  .;.'                     #
#               (,__....----'''       (,..--''                        #
#                                                                     #
#######################################################################

PCAT (Principal Component Analysis for Transients) is a suite of Python
utilities built to identify, isolate and characterize and classify 
transients found in LIGO/Advanced LIGO subsystems.

PCAT is built on standard Python libraries, such as numpy, scipy, matplotlib.
The scikit-learn module (http://scikit-learn.org/stable/, 
http://pypi.python.org/pypi/scikit-learn/) is also required. gwpy
(https://github.com/gwpy) is  also required when running PCAT_configread to
retrieve locked segments. pylal is required for data retrieval (download_frames.py)

*******************************************************************************
***************************      Installation:      ***************************
*******************************************************************************

To get started, add the PCAT folder to your PATH and PYTHONPATH
(or move PCAT.py and the necessary files in a folder in your PATH, PYTHONPATH).

To add the PCAT folder to your path simply add the following lines in your
.bashrc or .bash_profile folder in your home directory (~)

    EXPORT PATH=$PATH:/path/to/PCAT/
    EXPORT PYTHONPATH=$PYTHONPATH:/path/to/PCAT/
    eval `/ligotools/bin/use_ligotools`
    
Either source scikit-learn from an existing location:

or install scikit-learn from scratch:
    http://scikit-learn.org/stable/, 
    http://pypi.python.org/pypi/scikit-learn/

gwpy also has to be sourced:
    source ~detchar/opt/gwpysoft/etc/gwpy-user-env.sh 
*******************************************************************************
***************************        Contents:        ***************************
*******************************************************************************

- Main Program, which runs all the pipeline and returns URL for the results:
    - PCAT.py
    - PCAT_configread.py  (batch submit script + summary pages)
  If running time-domain analysis, the given time interval is scanned for 
  transients, on which PCA is performed, redusing these to a set of Principal
  Component scores, are then clustered.
  If running frequency-domain analysis, the given time interval is split into
  subsegments of which PSDs are computed. These are again clustered in the
  principal components space and divided by type.
  
  Results consist of:
   - Clickable scatterplots of the first few principal component scores (by
     default the first three) for each observation;
   - Clickable glitchgrams of the analyzed time interval:
   - Triangle plot containing all the first 15 components scores scatterplots
     plus single component distributions.
   - Breakdown of observations by type;
   - Representation of the average observation of each type in time/frequency
     domain (depending on the type of analysis being performed);
   - Database (pickled python list) containing all the observations of the run
     plus useful metadata, such as GPS peak, SNR and peak frequency for time
     domain and GPS start/end for the segment in frequency domain analysis;
   - List of the segments analyzed and total analyzed time.
  Also, if doing a time-domain analysis, results include:
   - Time series of identified transients, showing both raw time series
     and principal component-reconstructed waveforms;
   - Spectrograms for the representative transients;
   - Glitchgram showing the time distribution of the transients.
  If doing frequency domain analysis, results include:
   - PSDs plots and principal component-reconstructed PSDs.
  
  Extra:
   - Plots of the first few (by default 10) principal components;
   - Matrix of the Principal Components, as a binary file (use with
     numpy.load());
   - Explained variance as a function of the number
     of Principal Components used in reconstructing the data.


All of the following contain definitions/functions used in PCAT.py, and
are *NECESSARY* for PCAT.py to work correctly, but they can also be run
standalone to run each step of the pipeline:
    - download_frames.py
    - data_conditioning.py
    - finder.py
    - PCA.py
    - GMM.py

The following contain *NECESSARY* functions/definitions used by the above:
     - spike_class.py
     - utilities_PCAT.py

- Extra utilities (require all all of the above to run):
    Data conditioning:
       - pickled_butter.py
       - PSD.py	
    Plotting:
       - spikes_time_series.py
       - time_series.py
    Data manipulation:
       - database.py

All of these (apart from 'spike_class.py' and 'utilities_PCAT.py') can be run 
from the command line:
	$ PCAT.py [-h] [arguments]
	$ time_series.py [-h] [arguments]	
or
	$ python PCAT.py [-h] [arguments]
	$ python time_series.py [-h] [arguments]

Get usage on each program by calling it with '-h' or '--help':
	
	$ PCAT.py -h
	$ GMM.py -h

*******************************************************************************
***********             Brief description of contents:              ***********
*******************************************************************************

All of the programs in PCAT have a built-in help describing usage.
To show usage: either run the program with no arguments or  '-h' or '--help'.

To quickly get started after installation, scroll down to the PCAT.py examples.
PCAT.py.

The PCAT pipeline:
    1) Retrieve and condition Data: 
        Retrieve time series from frame files using lalframe, then
        prepare data for analysis rurnning through the conditioning
        functions, one in the following:
        whiten(...), butterworth_band_pass(), compute_psd(...)
        (from data_conditioning.py)
        Before actually conditioning data, PCAT looks for existing
        files and tries to load them, in order to speed up pipeline
        if already run before

    3) Create database from conditioned files:
        create_data_matrix(...), create_data_matrix_from_psd(...)
        from utilities_PCAT.py
    
    4) Perform PCA, get scores matrix and Principal Components:
        PCA(...) from PCA.py
        
    5) Cluster the scores matrix from PCA.py using a small number of
       Principal Components (user-selectable):
       gaussian_mixture(...) from GMM.py
        
       Another clustering algorithm can be implemented by changing a
       few lines of code (once one has defined the clustering algorithm )
       defined. See GMM.py
    
    6) Plot clickable scatterplots, clickable glitchgram, time series or PSDs:
        scatterplot(...), plot_glitchgram(), plot_time_series(...), plot_psds(...)
        from GMM.py and utilities_PCAT.py
    
    7) Print URL to analysis results.
        Output file (in the output directory) is a ".list" binary file.
        This is simply a python list of Spike() objects, which can be
        loaded into numpy using pickle's load() function (spike_class
        has to be imported for this to work).

PCAT.py's output is an URL with scatterplots and analysis for the given times.


    - PCAT.py               Full pipeline for a single channel.
    - PCAT_configread.py    Wraps the above and runs on a list of channels,
                            generating html summary pages.
                            This can take as argument either start and end GPS times,
                            resulting in the analysis of the locked times between
                            start and end times or a list of GPS times to be analyzed.
                            Configuration files are also needed (see .config files
                            for examples).
    - download_frames.py    Retrieves frame files using gwpy library
                            Arguments are start and end time in GPS time, 
                            channel name, IFO and frame type.
                            Data is downloaded and split into
                            segments for easier manipulation.
                            Default segment size is 60 seconds.
                            *** Does NOT download data if files already exist***
                             
    - PCA.py                Perform PCA on input data and plots the scores
                            for the first few components.
                            Also plots the explained variance vs. the number
                            of Principal Components. The number of plotted
                            components can be chosen through '--extras'.
                            
                            Argument(s): a single (or list of) output(s) from
                            finder.py or database.py.
                            
    - GMM.py                Performs PCA on input data and clusters in the
                            Principal Component space.
                            Input is a single (or list of) output(s) from
                            finder.py or database.py
                            GMM can be called with the '-r' option to remove
                            certain clusters after PCA and clustering have been
                            completed. 
                            A new database file without the observations
                            corresponding to the removed clusters is saved.
                            GMM.py can be then run again on the new database.
                            This is useful to remove unwanted clusters and/or                             
                            observations from a database.
                            
    - finder.py             Inputs either pickled or plain text and searches
                            for transients. Output is a list of instances of the
                            Spike class containing details about the transient,
                            (see spike_class.py).
                            Each sampled transient is centered on the maximum of
                            the absolute value of the points.
                            A simple way to determine the threshold is to use
                            
                               time_series.py -t threshold file1 file2 file3 ..
                            
                            where file1 file2 file3 ... are the names of the
                            files containing the time series.
                            to be analyzed and 'threshold' is a float
                            corresponding to trigger threshold in units of the
                            standard deviation of each file.
                            Two threshold lines are plotted at 
                            +threshold*sigma and -threshold*sigma.
                            
Extra utilities:
     - database.py        Create database files readable by PCA.py and
                          GMM.py from either plain-text or pickled files.
                          This can be also used to merge already
                          existing databases, see -h.
                          
     - pickled_butter.py  A fourth order back and forth butterworth
                          band-pass filter, output is a pickled file which 
                          can be used in GMM.py and PCA.py.
                          Argument: **One** plain text file. 
                          (use: with xargs -P  )
                          
     - time_series.py     Plots a time series the input files and two threshold
                          lines, defined by the -t N parameter, corresponding to
                          +N*sigma and -N*sigma, threshold in units of the
                          standard deviation of the input file.
                          ** WARNING: this discards 30% of the
                          time series from the beginning and the end, thus only
                          showing 40% of the length of the actual time series
                          This is used to avoid plotting ringing artifacts due
                          to filtering when filtering files.
                          Use --remove_seconds 0 to avoid this (for example
                          when plotting raw data)
                          
     - spikes_time_series.py Plots time series of the transients in a database
                             output from finder.py.
                          
     - PSD.py                Takes as input single or multiple pickled or plain
                             text files and returns a pickled numpy array
                             containing the PSD.
                             Plot the PSD using --plot.
Misc:
     - spike_class.py contains Spike() class definitions.
                     A Spike() object is used to store information about the
                     transients found with finder.py and is also used used in 
                     GMM.py, PCA.py, merge_database.py and
                     spikes_time_series.py.
                     Information stored in a Spike() instance:
                        spike.start         -> Transient starting point.
                        spike.end           -> Transient ending point
                        spike.width         -> Number of points above N*sigma.
                                               (spike.end-spike.start)
                        spike.peak          -> Transient peak point.
                        spike.norm          -> Transient normalization constant
                                               (for use with PCA.py).
                        spike.waveform      -> Time Series
                        spike.segment_start -> Segment start expressed 
                                               in GPS time.
                        spike.segment_end   -> Segment end expressed
                                               in GPS time.
                        spike.peak_GPS      -> Spike peak expressed 
                                               in GPS time.
                        spike.SNR           -> SNR of the glitch (time domain)
                                               This is only meaningful when data
                                               is whitened.
                                               
                                             
     - utilities_PCAT.py Contains various functions (and module imports) used
                         all throughout PCAT.


*******************************************************************************
*****************          Usage (quick start)               ******************
*******************************************************************************

Quickly get started:
PCAT can be used to find noise transients in the time domain and classify and
characterize them, or it can be used to classify and characterize power spectral
densities of short segments.

Each time series is read from frame files and, if requested (off by default),
is saved (in a binary format readable using numpy's load() function) to:
    ~/PCAT/Data/${CHANNEL_NAME}/${INTERVAL_IDENTIFIER}/raw_data
where ${INTERVAL_IDENTIFIER} is either the name of the supplied list or start
time and end time as supplied from the command line.

The time series is processed and saved to
    ~/PCAT/Data/${CHANNEL_NAME}/${INTERVAL_IDENTIFIER}/${CONDITIONED_FOLDER}
where again where again ${CONDITIONED_FOLDER} depends on the command line
arguments and the files are again binary, readable with numpy's load().

A database of the is created and saved to the output folder, either:
~/public_html/time_PCAT/${CHANNEL_NAME}/${INTERVAL_IDENTIFIER}/${PARAMETERS}
or
~/public_html/frequency_PCAT/${CHANNEL_NAME}/${INTERVAL_IDENTIFIER}/${PARAMETERS}

An URL to the folder where the database and plots are saved is returned.

PCAT can also be run on a batch of channels using PCAT_configread.py, which takes
as arguments either a list of times or a couple GPS times and two configuration
files, one for time domain and one for frequency domain (see .config files for
examples) and returns a summary page containing links to results for each channel.

See PCAT_configread.py -h for more.

Usage:

- Time Domain Analysis:
     PCAT.py --size 8 --padding_seconds 4 --time --start 965368792 \
                --end 965368999 --whiten -v 500 --IFO L \
                -c L1:LSC-DARM_ERR --frame R -t 5.5
     
     Runs PCAT, time-domain analysis. PCA is performed on transients sampled
     using a number of points given by -v option. Small values (less than 8000)
     gives faster results. The temporal length of these transients in seconds is
     number_points/sampling_frequency. ** WARNING: '--whiten' downsamples input
     data ** to 4096 Hz, so in this case we have 500/4096=0.122s per transient.
     
     In this case the analysis is performed from GPS times 965368792 to
     965368999, analyzed in 8 seconds chunks (downloading 8 seconds segments
     and padding with 4 seconds at beginning the start of the chunks, in order
     to avoid ringing artifacts due to filtering and FFT) and then applying a
     whitening filter.
     Channel name is L1:LSC-DARM_ERR with sampling frequency 16384 Hz,
     from Livingston Observatory (--IFO L), frame type 'R' (frame type can be
     found using ligo_data_find).
     
     The --whiten whitens the data before running the trigger finder.
     Data is also downsampled to 4096 Hz (ANALYSIS_FREQUENCY, defined in
     data_conditioning.py) before whitening for faster processing.
     
     Trigger threshold is set at 5.5 sigma through -t. This can be found by
     trial and error, or simply plotting some of the conditioned time series
     using time_series.py with the -t option and looking at the plots for
     different values of the threshold.

     
     The --size 8 option sets the length of the analyzed chunks to 8
     seconds. These chunks are padded with 4 seconds at beginning and the end,
     which are then discarded. 
     This is done to exclude ringing artifacts due to the Discrete Fourier
     Transforms and the edges of the segments. 
     
     The shorter --size segment_size the faster the analysis (whitening and
     filtering take long for longer segments).
     
     When changing --size one should also use --padding_overlap to change the 
     percentage by which segments are overlapped and consequently, the amount
     of time excluded from the beginning and end of the segments. 
     
- Frequency Domain Analysis:
     PCAT.py --frequency -v 8192 --start 965282415 --end 965282915 --frame R
                    -I L -c L1:OMC-PD_SUM_OUT_DAQ
     
     PCA is performed in the frequency domain, using power spectral densities as
     observations.
     
     The only difference with the above example is the --frequency flag, and the
     absence of the filter flags (--filter, --low, --high).
     
     This time -v sets the resolution of the power spectral density. For faster
     processing, a power of 2 should be chosen (e.g. 4096, 8192, 16384, ...),
     the resolution of the PSD is given by nyquist_frequency/variable_number,
     where the nyquist frequency is half of the sampling frequency
     
     In this case the sampling frequency is 32768, and variable number is set to
     8192, yielding a PSD resolution of 4 Hz.
     
     One can use the --low low_frequency and --high high_frequency options to
     compute a 1 Hz PSD and perform PCA using as observations only the band
     between low_frequency and high_frequency of the power spectrum.
     In this case the number of variables used is high_frequency-low_frequency.
     This means that using '--low 30 --high 2000' sets -v to 1970.
     Keep in mind that the higher -v, the slower PCA decomposition is.


*******************************************************************************
*************         Usage (output database)                  ****************
*******************************************************************************

PCAT outputs pickled list of Spike() istances. These can be used with
PCA() and GMM() to avoid running the full pipeline twice:
Usage example:
        PCA.py --time t-5.5_w-1500.list
or
        GMM.py -s 32768 --time t-5.5_w-1500.list
  
For frequency-domain:

        PCA.py --frequency psd_database.data
or
        GMM.py -s 32768 --frequency psd_database.data

For usage GMM.py -h and PCA.py -h.


Output files can be merged using database.py (call with -h for help).

