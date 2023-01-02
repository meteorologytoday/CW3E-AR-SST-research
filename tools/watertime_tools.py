from datetime import datetime, timedelta
import numpy as np



def getWateryear(d):

    return d.year + 1 if d.month in [10, 11, 12] else d.year


def getWatertime(d):

    wateryear = getWateryear(d)
    beg_of_wateryear = datetime(wateryear-1, 10, 1)
    end_of_wateryear = datetime(wateryear, 4, 1)
   
    normalized_time = (d - beg_of_wateryear) / (end_of_wateryear - beg_of_wateryear)

    return normalized_time
    
