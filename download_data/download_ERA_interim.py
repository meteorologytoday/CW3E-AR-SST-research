#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1979-01-01/to/1979-01-02",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "34.128/165.128/166.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
    "target": "output",
})

server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "1979-01-01/to/1979-01-02",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "34.128/165.128/166.128",
    "step": "3/6/9/12",
    "stream": "oper",
    "time": "00:00:00/12:00:00",
    "type": "fc",
    "target": "output",
})
