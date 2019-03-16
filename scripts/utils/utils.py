"""
Set of utility functions
"""

#http://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
def get_overlap(a,b):
    """
    Compute the overlap between two discrete and closed intervals (borders are inclusive) a and b
    :param a: An interval written as a list. E.g. [5,15]
    :param b: An interval written as a list. E.g. [10,20]
    :return: The size of the overlap
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)