#!/usr/bin/env python
import sys
import os
from argparse import ArgumentParser

import numpy as np

import lal
from glue.ligolw import ligolw, utils, lsctables
import cPickle as pickle


argp = ArgumentParser(description="Convert PCAT information into an appropriate set of LIGO_LW tables.")
argp.add_argument('--database', help="File to process. Required, (PCAT '.list.' file)")
argp.add_argument('--channel', help="Name of channel, required. E.g. 'H1:FAKE-STRAIN'")
argp.add_argument('--search', default='PCAT', help="Name of search, default is 'PCAT'.")
args = argp.parse_args()

try:
    ifo, channel = args.channel.split(":")
except ValueError:
    sys.exit("Invalid channel specification")

proc_table = lsctables.New(lsctables.ProcessTable, ['comment', 'process_id', 'ifos', 'username', 'program'])
proc_row = proc_table.RowType()
proc_row.process_id = pid = proc_table.get_next_id()
proc_row.ifos = set([ifo])
proc_row.program = sys.argv[0]
proc_row.username = os.path.expanduser("~").split("/")[-1]
proc_row.comment = "converted from PCAT output"
proc_table.append(proc_row)

#TODO: ADD PCAT RESULTS URL HERE

table = lsctables.New(lsctables.SnglBurstTable, ["ifo", "peak_time", 
      "peak_time_ns", "start_time", "start_time_ns",
      "duration",  "search", "event_id", "process_id",
      "central_freq", "channel", "amplitude", "snr",
      "bandwidth"]) #"confidence","chisq", "chisq_dof"

#TODO add PCAT type to table??

with open(args.database, "rb") as f:
    database = pickle.load(f)

class GPS():
    def __init__(self, seconds, decimals):
        self.seconds = int(seconds)
        self.nanoseconds = int(decimals*1e8)
    
for index, spike in enumerate(database):
    row = table.RowType()
    
    tmp_time = int(np.floor(spike.peak_GPS))
    peak_GPS = GPS(tmp_time, spike.peak_GPS-tmp_time)
    row.set_peak(peak_GPS)
    
    peak_tmp = spike.peak_GPS + spike.start*spike.sampling
    peak_tmp_int = int(np.floor(peak_tmp))
    start = GPS(peak_tmp_int, peak_tmp-peak_tmp_int)
    row.set_start(start)
    # duration is temporary, replace with spike.duration
    row.duration = (spike.waveform.size)/float(spike.sampling)
    #row.duration = spike.duration
    
    # Duration is simply f_max -f_min, or nyquist frequency minus lowest frequcency (1/segment_len)
    row.bandwidth = spike.sampling/2.0 - spike.sampling/float(spike.waveform.size)
    #FIXME: Bandwidth is actually lower because data is high pass filtered.
    # We cannot do much about this until we add this to PCAT.py
    
    ## Confidence is simply PCAT's threshold (LIGO-T1200125), either add this to PCAT database (add in finder.py)
    ## or simply wait to merge xml_postproc into PCAT as the above
    #row.confidence = 0.0 # FIXME: Do we care about this right?
    row.central_freq = spike.central_freq

    #row.chisq = 0
    #row.chisq_dof = 2*band*dur
    

    row.snr = spike.SNR
    row.amplitude = np.max(np.abs(spike.waveform))
    
    #row.type =  spike.type
    row.ifo = ifo
    row.channel = channel
    row.search = args.search
    
    row.event_id = table.get_next_id()
    row.process_id = pid

    table.append(row)

# Build output document
xmldoc = ligolw.Document()
xmldoc.appendChild(ligolw.LIGO_LW())
xmldoc.childNodes[0].appendChild(proc_table)
xmldoc.childNodes[0].appendChild(table)

out_name = args.database.replace('.list', '.xml.gz')
utils.write_filename(xmldoc, out_name, gz=True)

print "Saved '{0}'".format(out_name)


# Done!