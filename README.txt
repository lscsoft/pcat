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
http://pypi.python.org/pypi/scikit-learn/) is also required.
gwpy (https://github.com/gwpy) is  also required when running pcat-multi to
retrieve locked segments. pylal is required for data retrieval (data.py)

*******************************************************************************
***************************      Installation:      ***************************
*******************************************************************************

To get started install PCAT with

pip install git+https://github.com/lscsoft/pcat.git

or for a local install

git clone https://github.com/lscsoft/pcat
cd pcat
pip install .

Then add the pcat/bin folder to your PATH and the pcat/ and pcat/pcat folders to your
PYTHONPATH in the .bashrc or .bash_profile files 

    EXPORT PATH=$PATH:/path/to/pcat/bin
    EXPORT PYTHONPATH=$PYTHONPATH:/path/to/pcat/:/path/to/pcat/pcat/

Source the Lab-LDG GWpy activate script:

. ~detchar/opt/gwpysoft/bin/activate

Local install  matplotlib and scipy

pip install --user --upgrade matplotlib
pip install --user --upgrade scipy

Then local install scikit-learn 

git clone https://github.com/scikit-learn/scikit-learn
cd scikit-learn
python setup.py install --user 

Then add the local sklearn folder to your PYTHONPATH in the .bashrc or .bash_profile files 

export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python2.7/site-packages/sklearn/

You should be ready to go...

*******************************************************************************
***************************        Contents:        ***************************
*******************************************************************************

- Main Program, which runs all the pipeline and returns URL for the results:
    - pcat
    - pcat-multi  (batch submit script + summary pages)
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


All of the following contain definitions/functions used in pcat, and
are *NECESSARY* for pcat to work correctly, but they can also be run
standalone to run each step of the pipeline:
    - data.py
    - condition.py
    - finder.py
    - pca.py
    - gmm.py

The following contain *NECESSARY* functions/definitions used by the above:
     - spike.py
     - utils.py

- Extra utilities (require all of the above to run):
    Data conditioning:
       - pickled_butter.py
       - PSD.py	
    Plotting:
       - spikes_time_series.py
    Data manipulation:
       - database.py

All of these (apart from 'spike.py' and 'utils.py') can be run 
from the command line:
	$ pcat [-h] [arguments]
	$ spikes_time_series.py [-h] [arguments]	
or
	$ python pcat [-h] [arguments]
	$ python spike_time_series.py [-h] [arguments]

Get usage on each program by calling it with '-h' or '--help':
	
	$ pcat -h
	$ pcat-multi -h
	$ gmm -h

*******************************************************************************
***********             Brief description of contents:              ***********
*******************************************************************************

All of the programs in pcat have a built-in help describing usage.
To show usage: either run the program with no arguments or  '-h' or '--help'.

To quickly get started after installation, scroll down to the pcat examples.

The pcat pipeline:
    1) Retrieves and condition Data: 
        Retrieves time series from frame files using lalframe, then
        prepares data for analysis rurnning through the conditioning
        functions, one in the following:
        whiten(...), butterworth_band_pass(), compute_psd(...)
        (from condition.py)
        Before actually conditioning data, pcat looks for existing
        files and tries to load them, in order to speed up pipeline
        if already run before

    3) Creates database from conditioned files:
        create_data_matrix(...), create_data_matrix_from_psd(...)
        from utils.py
    
    4) Performs the PCA, get scores matrix and Principal Components:
        PCA(...) from pca.py
        
    5) Clusters the scores matrix from pca.py using a small number of
       Principal Components (user-selectable):
       gaussian_mixture(...) from gmm.py
        
       Another clustering algorithm can be implemented by changing a
       few lines of code (once one has defined the clustering algorithm )
       defined. See gmm.py
    
    6) Plots clickable scatterplots, clickable glitchgram, time series or PSDs:
        scatterplot(...), plot_glitchgram(), plot_time_series(...), plot_psds(...)
        from pcat.gmm and pcat.utils
    
    7) Prints URL to analysis results.
        Output file (in the output directory) is a ".list" binary file.
        This is simply a python list of Spike() objects, which can be
        loaded into numpy using pickle's load() function (spike_class
        has to be imported for this to work).

pcat's output is an URL with scatterplots and analysis for the given times.

    - pcat                  Full pipeline for a single channel.
    - pcat-multi            Wraps the above and runs on a list of channels,
                            generating html summary pages.
                            This can take as argument either start and end GPS times,
                            resulting in the analysis of the locked times between
                            start and end times or a list of GPS times to be analyzed.
                            Configuration files are also needed (see .config files
                            for examples).
    - data.py               Retrieves frame files using gwpy library
                            Arguments are start and end time in GPS time, 
                            channel name, IFO and frame type.
                            Data is downloaded and split into
                            segments for easier manipulation.
                            Default segment size is 60 seconds.
                            *** Does NOT download data if files already exist***
                             
    - pca.py                Performs PCA on input data and plots the scores
                            for the first few components.
                            Also plots the explained variance vs. the number
                            of Principal Components. The number of plotted
                            components can be chosen through '--extras'.
                            
                            Argument(s): a single (or list of) output(s) from
                            finder.py or database.py.
                            
    - gmm.py                Performs PCA on input data and clusters in the
                            Principal Component space.
                            Input is a single (or list of) output(s) from
                            finder.py or database.py
                            gmm can be called with the '-r' option to remove
                            certain clusters after PCA and clustering have been
                            completed. 
                            A new database file without the observations
                            corresponding to the removed clusters is saved.
                            gmm.py can be then run again on the new database.
                            This is useful to remove unwanted clusters and/or                             
                            observations from a database.
                            
    - finder.py             Inputs either pickled or plain text and searches
                            for transients. Output is a list of instances of the
                            Spike class containing details about the transient,
                            (see spike.py).
                            Each sampled transient is centered on the maximum of
                            the absolute value of the points.
                            A simple way to determine the threshold is to use
                            
                               spikes_time_series.py -t threshold file1 file2 file3 ..
                            
                            where file1 file2 file3 ... are the names of the
                            files containing the time series.
                            to be analyzed and 'threshold' is a float
                            corresponding to trigger threshold in units of the
                            standard deviation of each file.
                            Two threshold lines are plotted at 
                            +threshold*sigma and -threshold*sigma.
                            
Extra utilities:
     - database.py        Create database files readable by PCA.py and
                          gmm.py from either plain-text or pickled files.
                          This can be also used to merge already
                          existing databases, see -h.
                          
     - pickled_butter.py  A fourth order back and forth butterworth
                          band-pass filter, output is a pickled file which 
                          can be used in gmm.py and pca.py.
                          Argument: **One** plain text file. 
                          (use: with xargs -P  )
                          
     - spikes_time_series.py Plots time series of the transients in a database
                             output from finder.py.
                          
     - PSD.py                Takes as input single or multiple pickled or plain
                             text files and returns a pickled numpy array
                             containing the PSD.
                             Plot the PSD using --plot.
Misc:
     - spike.py contains Spike() class definitions.
                     A Spike() object is used to store information about the
                     transients found with finder.py and is also used used in 
                     gmm.py, pca.py, database.py and
                     spikes_time_series.py.
                     Information stored in a Spike() instance for time domain:
                        spike.start             -> Transient starting point.
                        spike.end               -> Transient ending point
                        spike.width             -> Number of points above N*sigma.
                                                  (spike.end-spike.start)
                        spike.peak              -> Transient peak point.
                        spike.norm              -> Transient normalization constant
                                                  (for use with PCA.py).
                        spike.waveform          -> Time Series
                        spike.segment_start     -> Segment start expressed 
                                                  in GPS time.
                        spike.segment_end       -> Segment end expressed
                                                  in GPS time.
                        spike.peak_GPS          -> Spike peak expressed 
                                                  in GPS time.
                        spike.SNR               -> SNR of the glitch (time domain)
                                                  This is only meaningful when data
                                                  is whitened.
                        spike.sampling          -> Sampling Frequency
                        spike.type              -> Type
                        spike.peak_frequency    -> Peak frequency  
                    For frequency domain:
                        spike.start             -> Segment start (GPS time)
                        spike.end               -> Segment end (GPS time)
                        spike.waveform          -> PSD for the segment
                        spike.sampling          -> Sampling frequency
                    
                    To get a list of all the defined attributes (python code): 
                    > import inspect
                    > variables = [i for i in dir(t) if not inspect.ismethod(i)
                
     - utils.py          Contains various functions (and module imports) used
                         all throughout PCAT.


*******************************************************************************
*****************          Usage (quick start)               ******************
*******************************************************************************

Quickly get started:
pcat can be used to find noise transients in the time domain and classify and
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

pcat can also be run on a batch of channels using pcat-multi, which takes
as arguments either a list of times or a couple GPS times and two configuration
files, one for time domain and one for frequency domain (see .config files for
examples) and returns a summary page containing links to results for each channel.

See pcat-multi -h for more.

Example of time-domain analysis:
     pcat  --silent --time --whiten --highpass 10 --IFO L \
           --frame L1_HOFT_C00 -c L1:GDS-CALIB_STRAIN -t 4.5 -v 1024 \
           --resample 8192 --components 20 --reconstruct \
           --start 1126224017 --end  1126224017
     Runs PCAT, time-domain analysis. PCA is performed on transients sampled
     using a number of points given by -v option. Small values (less than 8000)
     gives faster results. The temporal length of these transients in seconds is
     number_points/sampling_frequency. 
     
     In this case the analysis is performed from GPS times 1126224017 to
     1126224017 with a highpass and a whitening filter. The channel is
     L1:GDS-CALIB_STRAIN with sampling frequency 16384 Hz, downsampled at 8192 Hz
     from Livingston Observatory (--IFO L), frame type 'L1_HOFT_C00'. Trigger threshold
     is set at 4.5 sigma through -t.
     
Example of frequency-domain analysis:
     pcat --frequency -v 8192 --start 965282415 --end 965282915 --frame R
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
        pca.py --time t-5.5_w-1500.list
or
        gmm.py -s 32768 --time t-5.5_w-1500.list
  
For frequency-domain:

        pca.py --frequency psd_database.data
or
        gmm.py -s 32768 --frequency psd_database.data

For usage gmm.py -h and pca.py -h.

This can be also loaded in python using the following python code:

> import pickle
> with open("database.list", "rb") as f:
>   database = pickle.load(f)

if the database is from time domain analysis,


".list" files can be merged using database.py (call with -h for help).

