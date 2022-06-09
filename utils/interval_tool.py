import numpy as np


class IntervalTree():
    def __init__(self,length=100000):
        self._interval = np.zeros((length, 2), int)
        self._current_index = 0

    def __repr__(self):
        return str(self._interval)

    def insert_interval(self, start, end):
        self._interval[self._current_index] = [start, end]
        self._current_index+=1

    def find_index(self, start, end):
        return np.argwhere((start < self._interval[:, 1]) & (self._interval[:, 0] < end))

    def find(self, start, end):
        return self._interval[self.find_index(start, end)[:, 0]]


if __name__ == "__main__":
    import random

    import log
    logger = log.get_logger("IntervalTree Test", log.DEBUG)
    intervals = IntervalTree()
    for i in range(10):
        start = random.randint(0, 100)
        end = start + random.randint(0, 10)
        intervals.insert_interval(start, end)
    query = (10, 50)
    logger.debug("query: {} {}".format(query[0], query[1]))
    logger.debug("intervals: {}".format(intervals))
    logger.debug("new intervals: {}".format(
        intervals.find(query[0], query[1])))
