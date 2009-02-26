__author__="fccoelho"
__date__ ="$26/02/2009 10:44:29$"

import Gnuplot
from numpy import histogram

class RTplot:
    '''
    Real time plotting class based on Gnuplot
    '''
    def __init__(self):
        self.gp = Gnuplot.Gnuplot(persist = 1)
    def plotlines(self,data,names=[]):
        '''
        Create a sinlge/multiple line plot from a numpy array or record array.
        :Parameters:
            - `data`: must be a numpy array or record array, with series as rows
            - `names`: is a list of strings to serve as legend labels
        '''
        self.gp('set data style lines')
        if isinstance(data,numpy.core.records.recarray):
            return self._linesFromRA(data)
        if len(data.shape) > 1 and len(data.shape) <= 2:
            plots = []
            for n,row in enumerate(data):
                plots.append(Gnuplot.PlotItems.Data(enumerate(row),title=names[n]))
            self.gp.plot(*tuple(plots))
        elif len(data.shape) >2:
                pass
        else:
            self.gp.plot(Gnuplot.PlotItems.Data(data,names[0]))

    def _linesFromRA(self,data):
        '''
        Record-array specific code
        '''
        plots = []
        for n in data.dtype.names:
            if len(data.shape) > 1 and len(data.shape) <= 2:
                i = 0
                for row in data[n]:
                    plots.append(Gnuplot.PlotItems.Data(data[n],title=n+':%s'%i))
                    i += 1
            elif len(data.shape) >2:
                pass
            # TODO: figure out what to do with higher dimensional data
            else:
                plots.append(Gnuplot.PlotItems.Data(data[n],title=n))
        self.gp.plot(*tuple(plots))

    def plothist(self,data, names=[]):
        '''
        Create a sinlge/multiple Histogram plot from a numpy array or record array.
        :Parameters:
            - `data`: must be a numpy array or record array, with series as rows
            - `names`: is a list of strings to serve as legend labels
        '''
        self.gp('set data style boxes')
        if isinstance(data,numpy.core.records.recarray):
            return self._histFromRA(data)
        if len(data.shape) > 1 and len(data.shape) <= 2:
            plots = []
            for n,row in enumerate(data):
                m,bins = histogram(row,normed=True,bins=50,new=True)
                d = zip(bins[:-1],m)
                plots.append(Gnuplot.PlotItems.Data(d,title=names[n]))
            self.gp.plot(*tuple(plots))
        elif len(data.shape) >2:
            pass
        else:
            m,bins = histogram(data,normed=True,bins=50,new=True)
            d = zip(bins[:-1],m)
            self.gp.plot(Gnuplot.PlotItems.Data(d,title=names[0]))

    def _histFromRA(self,data):
        '''
        Record-array specific code
        '''
        plots = []
        for n in data.dtype.names:
            if len(data.shape) > 1 and len(data.shape) <= 2:
                i = 0
                for row in data[n]:
                    m,bins = histogram(row,normed=True,bins=50,new=True)
                    d = zip(bins[:-1],m)
                    plots.append(Gnuplot.PlotItems.Data(d,title=n+':%s'%i))
                    i += 1
            elif len(data.shape) > 2:
                pass
            else:
                m,bins = histogram(data[n],normed=True,bins=50,new=True)
                d = zip(bins[:-1],m)
                plots.append(Gnuplot.PlotItems.Data(d,title=n))
        self.gp.plot(*tuple(plots))

if __name__ == "__main__":
    print "Hello";