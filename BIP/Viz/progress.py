__author__ = 'fccoelho'
import sys

def bar(count, total, suffix=""):
    """
    Prints a progress bar in standard out
    :param count: steps taken
    :param total: total steps
    :param suffix: Comment
    :return:
    """
    if suffix == "":
        suffix = 'Completed {} out of {}'.format(total, count)
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))