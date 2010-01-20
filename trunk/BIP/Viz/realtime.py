__author__="fccoelho@gmail.com"
__date__ ="$26/02/2009 10:44:29$"
__docformat__ = "restructuredtext en"
import Gnuplot
import numpy
import matplotlib
import pylab as P


class RTplot:
    '''
    Real time plotting class based on Gnuplot
    '''
    def __init__(self, persist=1,debug=0):
        self.gp = Gnuplot.Gnuplot(persist = persist, debug=debug)
        self.plots = []

    def clearFig(self):
        '''
        Clears the figure.
        '''
        self.plots = []
        #self.gp.reset()

    def scatter(self,x,y,names=[],title='',style='points'):
        """
        Makes scatter plots from numpy arrays.
        if arrays are multidimensional, multiple scatter plots will be generated, pairing rows.
        """
        if isinstance(x, numpy.ndarray):
            if not isinstance(y, numpy.ndarray):
                raise TypeError("If x is a numpy array, y must also be an array.")
            if x.shape != y.shape:
                raise ValueError("x, %s and y, %s arrays must have the same shape."%(x.shape,y.shape))
            if names:
                if len(names) != x.shape[0]:
                    raise ValueError("names list must have exactly %s items, but has %s."%(x.shape[0],len(names)))
        else:
            x = numpy.array(x)
            y = numpy.array(y)


        self.gp('set title "%s"'%title)
        if not names:
            names = ['s%s'%i for i in range(x.shape[0])]
        if len(x.shape) > 1 and len(x.shape) <= 2:
            i = 0
            for n in range(x.shape[0]):
                self.plots.append(Gnuplot.PlotItems.Data(x[n],y[n],title=names[i],with_=style))
                i += 1
            self.gp.plot(*tuple(self.plots))
        elif len(x.shape) >2:
                pass
        else:
            #print data
            self.plots.append(Gnuplot.PlotItems.Data(x,y,title=names[0],with_=style))
            self.gp.plot(*tuple(self.plots))
        
    def plotlines(self,data,names=[],title='',style='lines'):
        '''
        Create a sinlge/multiple line plot from a numpy array or record array.
        
        :Parameters:
            - `data`: must be a numpy array or record array, with series as rows
            - `names`: is a list of strings to serve as legend labels
            - `style`: plot styles from gnuplot: lines, boxes, points, linespoints, etc.
        '''
        self.gp('set title "%s"'%title)
        if isinstance (data, list):
            data = numpy.array(data)
        if isinstance(data,numpy.core.records.recarray):
            return self._linesFromRA(data,style)
        if len(data.shape) > 1 and len(data.shape) <= 2:
            i = 0
            for row in data:
                if names:
                    self.plots.append(Gnuplot.PlotItems.Data(row,title=names[i],with_=style))
                else:
                    self.plots.append(Gnuplot.PlotItems.Data(row,with_=style))
                i += 1
            self.gp.plot(*tuple(self.plots))
        elif len(data.shape) >2:
                pass
        else:
            #print data
            self.plots.append(Gnuplot.PlotItems.Data(data,title=names[0],with_=style))
            self.gp.plot(*tuple(self.plots))

    def _linesFromRA(self,data, style):
        '''
        Record-array specific code
        '''
        for n in data.dtype.names:
            if len(data.shape) > 1 and len(data.shape) <= 2:
                i = 0
                for row in data[n]:
                    self.plots.append(Gnuplot.PlotItems.Data(row,title=n+':%s'%i,with_=style))
                    i += 1
            elif len(data.shape) >2:
                pass
            # TODO: figure out what to do with higher dimensional data
            else:
                self.plots.append(Gnuplot.PlotItems.Data(data[n],title=n,with_=style))
        self.gp.plot(*tuple(self.plots))

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
        if not names:
            names = ['s%s'%i for i in range(data.shape[0])]
        if len(data.shape) > 1 and len(data.shape) <= 2:
            for n,row in enumerate(data):
                m,bins = numpy.histogram(row,normed=True,bins=50)
                d = zip(bins[:-1],m)
                self.plots.append(Gnuplot.PlotItems.Data(d,title=names[n]))
            self.gp.plot(*tuple(self.plots))
        elif len(data.shape) >2:
            pass
        else:
            m,bins = numpy.histogram(data,normed=True,bins=50)
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
                    m,bins = numpy.histogram(row,normed=True,bins=50)
                    d = zip(bins[:-1],m)
                    self.plots.append(Gnuplot.PlotItems.Data(d,title=n+':%s'%i))
                    i += 1
            elif len(data.shape) > 2:
                pass
            else:
                m,bins = numpy.histogram(data[n],normed=True,bins=50)
                d = zip(bins[:-1],m)
                self.plots.append(Gnuplot.PlotItems.Data(d,title=n))
        self.gp.plot(*tuple(self.plots))

class RTpyplot:
    """
    Creates an animated plot based on pylab
    """
    def __init__(self,  nseries=1, leng =10, names=['l'], title = ""):
        self.nseries = nseries
        self.leng = leng
        self.names = names
        self.title = title
        self.lines = []
        self.started = False

    def _setup(self):
        self.ax = P.subplot(111)
        self.canvas = self.ax.figure.canvas
        for i in range(self.nseries):
            line, = P.plot(numpy.arange(self.leng), [1]*self.leng,label=self.names[i],  title=self.title, animated=True)
            self.lines.append(line)
        P.legend(loc=0)
        P.grid()
        P.show()
        self.started = True
    def plotlines(self, data):
        if not self.started:
            self._setup()
        nlines = data.shape[0]
        assert nlines == self.nseries
        for i, d in enumerate(data):
            background = self.canvas.copy_from_bbox(self.ax.bbox)
            self.lines[i].set_ydata(d)
            self.ax.draw_artist(self.lines[i])
            self.canvas.blit(self.ax.bbox)

        
        
        
    
if __name__ == "__main__":
    pass
