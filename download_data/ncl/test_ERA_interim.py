#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1981-03-01/to/1981-03-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "146.128/147.128/165.128/166.128/167.128/176.128/177.128/182.128/228.128",
    "step": "6",
    "stream": "oper",
    "time": "00:00:00/12:00:00",
    "type": "fc",
    "target": "output2",
    "format": "netcdf",
})
