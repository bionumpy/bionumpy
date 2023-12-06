import itertools

import numpy as np

import bionumpy as bnp
from bionumpy import replace, LocationEntry
from bionumpy.string_array import StringArray


def find_indices(locations, intervals: bnp.Interval):
    """Find the indices of the locations in the intervals

    Parameters
    ----------
    locations : bnp.LocationEntry
        Locations
    intervals : bnp.Interval
        Intervals

    Returns
    -------
    np.ndarray
        Indices
    """
    interval_start = np.searchsorted(locations, intervals.start, side='left')
    interval_stop = np.searchsorted(locations, intervals.stop, side='right')
    n_locations = interval_stop - interval_start
    interval_indices = np.repeat(np.arange(len(intervals)), n_locations)
    change_indices = np.insert(np.cumsum(n_locations)[:-1], 0, 0)
    a = np.arange(len(interval_indices))-np.repeat(change_indices-interval_start, n_locations)
    return a, interval_indices

def map_locations(locations: LocationEntry, intervals):
    """Map locations to intervals

    Parameters
    ----------
    locations : bnp.LocationEntry
        Locations
    intervals : bnp.Interval
        Intervals

    Returns
    -------
    bnp.Interval
        Mapped intervals
    """
    location_indices, interval_indices = find_indices(locations.position, intervals)
    new_entries = locations[location_indices]

    names = intervals.name if hasattr(intervals, 'name') else StringArray(np.arange(len(intervals)).astype('S'))
    return replace(new_entries, chromosome=names[interval_indices], position=new_entries.position-intervals.start[interval_indices])

'''

tverrfaglig modellering
* entomologi
* metreologi
* epidemiologi
'''

