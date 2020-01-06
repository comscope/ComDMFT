import os
import glob
import warnings


def get_wiencase():
    # Determine WIEN2k case name
    case = os.path.basename(os.getcwd())
    # directory name is not case (happens when submitting to cluster)
    if not os.path.isfile(case + '.struct'):
        files = glob.glob('*.struct')
        if len(files) < 1:
            raise Exception, 'No struct file present.'
        elif len(files) > 1:
            # heuristic algorithm to determine case:
            # the struct file whose head is the same as the most files in this
            # directory is most likely the case
            candidates = [os.path.splitext(f)[0] for f in files]
            allheads = [os.path.splitext(f)[0] for f in os.listdir('.')]
            counts = [len([1 for f in allheads if f == candidate])
                      for candidate in candidates]
            index = counts.index(max(counts))
            case = candidates[index]
            warnings.warn(
                " case heuristically selected from multiple possibilities!")
        else:
            case, ext = os.path.splitext(os.path.basename(files[0]))
    return case
