__author__="fccoelho@gmail.com"
__date__ ="$26/02/2009 10:44:29$"
__docformat__ = "restructuredtext en"
import Gnuplot
import numpy


class RTplot:
    '''
    Real time plotting class based on Gnuplot
    '''
    def __init__(self):
        self.gp = Gnuplot.Gnuplot(persist = 1)
        self.plots = []

    def clearFig(self):
        '''
        Clears the figure.
        '''
        self.plots = []
        self.gp.reset()
    def plotlines(self,data,names=[],title='',style='lines'):
        '''
        Create a sinlge/multiple line plot from a numpy array or record array.
        :Parameters:
            - `data`: must be a numpy array or record array, with series as rows
            - `names`: is a list of strings to serve as legend labels
            - `style`: plot styles from gnuplot: lines, boxes, points, linespoints, etc.
        '''
        self.gp('set title "%s"'%title)
        if isinstance(data,numpy.core.records.recarray):
            return self._linesFromRA(data,style)
        if len(data.shape) > 1 and len(data.shape) <= 2:
            i = 0
            for row in data:
                self.plots.append(Gnuplot.PlotItems.Data(row,title=names[i],with_=style))
                i += 1
            self.gp.plot(*tuple(self.plots))
        elif len(data.shape) >2:
                pass
        else:
            self.plots.append(Gnuplot.PlotItems.Data(data,names[0],with_=style))
            self.gp.plot(*tuple(self.plots))

    def _linesFromRA(self,data, style):
        '''
        Record-array specific code
        '''
        for n in data.dtype.names:
            if len(data.shape) > 1 and len(data.shape) <= 2:
                i = 0
                for row in data[n]:
                    plots.append(Gnuplot.PlotItems.Data(row,title=n+':%s'%i,with_=style))
                    i += 1
            elif len(data.shape) >2:
                pass
            # TODO: figure out what to do with higher dimensional data
            else:
                plots.append(Gnuplot.PlotItems.Data(data[n],title=n,with_=style))
        self.gp.plot(*tuple(plots))

    def plothist(self,data, title='', names=[]):
        '''
        Create a sinlge/multiple Histogram plot from a numpy array or record array.
        :Parameters:
            - `data`: must be a numpy array or record array, with series as rows
            - `names`: is a list of strings to serve as legend labels
        '''
        self.gp('set data style boxes')
        self.gp('set title "%s"'%title)
        if isinstance(data,numpy.core.records.recarray):
            return self._histFromRA(data)
        if len(data.shape) > 1 and len(data.shape) <= 2:
            for n,row in enumerate(data):
                m,bins = numpy.histogram(row,normed=True,bins=50,new=True)
                d = zip(bins[:-1],m)
                self.plots.append(Gnuplot.PlotItems.Data(d,title=names[n]))
            self.gp.plot(*tuple(self.plots))
        elif len(data.shape) >2:
            pass
        else:
            m,bins = numpy.histogram(data,normed=True,bins=50,new=True)
            d = zip(bins[:-1],m)
            self.plots.append(Gnuplot.PlotItems.Data(d,title=names[0]))
            self.gp.plot(*tuple(self.plots))

    def _histFromRA(self,data):
        '''
        Record-array specific code
        '''
        for n in data.dtype.names:
            if len(data.shape) > 1 and len(data.shape) <= 2:
                i = 0
                for row in data[n]:
                    m,bins = numpy.histogram(row,normed=True,bins=50,new=True)
                    d = zip(bins[:-1],m)
                    self.plots.append(Gnuplot.PlotItems.Data(d,title=n+':%s'%i))
                    i += 1
            elif len(data.shape) > 2:
                pass
            else:
                m,bins = numpy.histogram(data[n],normed=True,bins=50,new=True)
                d = zip(bins[:-1],m)
                self.plots.append(Gnuplot.PlotItems.Data(d,title=n))
        self.gp.plot(*tuple(plots))

if __name__ == "__main__":
    print "Hello";