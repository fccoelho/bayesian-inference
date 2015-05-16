__author__ = 'fccoelho'
import sys
import time
from collections import OrderedDict


class MultiProgressBars:
    def __init__(self, bar_length=60):
        self.bars = OrderedDict()
        self.bar_length = bar_length
        self.clear()

    def __call__(self, barname, count, total, comment=""):
        """
        Updates the bars
        :param barname: Name of the bar to update
        :param count:
        :param total:
        :param comment:
        :return:
        """
        self.clear()
        self.bars[barname] = (count, total, comment)

        for barn in self.bars:
            count, total, comment = self.bars[barn]
            if comment == "":
                Comment = 'Completed {} out of {}'.format(self.bars[barn][0], self.bars[barn][1])
            else:
                Comment = comment if barn == barname else 'Completed {} out of {}'.format(self.bars[barn][0],
                                                                                          self.bars[barn][1])
            filled_len = int(round(self.bar_length * count / float(total)))
            percents = round(100.0 * count / float(total), 1)
            bar = '=' * filled_len + '-' * (self.bar_length - filled_len)
            sys.stdout.write('{}:[{}] {}{} ...{}\n'.format(barn, bar, percents, '%', Comment))

    def clear(self):
        """Clear screen, return cursor to top left"""
        sys.stdout.write('\033[2J')
        sys.stdout.write('\033[H')
        sys.stdout.flush()

    def refresh(self):
        if len(self.bars) == 0:
            return self.clear()
        barname = list(self.bars.keys())[0]
        count, total, comment = self.bars[barname]
        self.__call__(barname, count, total, comment)

if __name__ == "__main__":
    pb = MultiProgressBars()
    for i, j in zip(range(100), range(0, 200, 5)):
        pb("first", i, 100, comment="{} of {}".format(i, 100))
        pb("second", j, 200)
        pb("third", i, 100)
        time.sleep(.2)
