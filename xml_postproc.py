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

table = lsctables.New(lsctables.SnglBurstTable, ["ifo", "peak_time", 
      "peak_time_ns", "start_time", "start_time_ns",
      "duration",  "search", "event_id", "process_id",
      "central_freq", "channel", "amplitude", "snr"])
      #"confidence","chisq", "chisq_dof", "bandwidth"])

#TODO add PCAT type to table??

with open(args.database, "rb") as f:
    database = pickle.load(f)
class GPS():
    def __init__(self, seconds, decimals):
        self.seconds = seconds
        self.nanoseconds = decimals*1e10
    
for index, spike in enumerate(database):
    row = table.RowType()
    
    tmp_time = str(spike.peak_GPS).split('.')
    peak_GPS = GPS(int(tmp_time[0]), int(tmp_time[1]))
    #row.set_peak(lal.LIGOTimeGPS((gps_seconds, gps_nanoseconds)))
    row.set_peak(peak_GPS)
    
    tmp_time = str(spike.peak_GPS+spike.start*spike.sampling).split('.')
    start = GPS(int(tmp_time[0]), int(tmp_time[1]))
    row.set_start(start)
    # duration is temporary, replace with spike.duration
    row.duration = (spike.end - spike.start)*spike.sampling 
    #row.duration = spike.duration
    
    #row.bandwidth = band
    
    row.central_freq = spike.central_freq

    #row.chisq = 0
    #row.chisq_dof = 2*band*dur
    #row.confidence = 0.0 # FIXME: Do we care to do this right?

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

with open(args.database.replace('.list', '.xml'), 'w') as output:
    xmldoc.write(fileobj=output)