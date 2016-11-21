from __future__ import absolute_import
from six.moves import range
from six.moves import zip
import numpy
import pylab as P
from six.moves.xmlrpc_server import SimpleXMLRPCServer
from multiprocessing import Process
from subprocess import Popen, PIPE

__author__ = "fccoelho@gmail.com"
__date__ = "$26/02/2009 10:44:29$"
__docformat__ = "restructuredtext en"

__ports_used = []


class RTplot:
    '''
    Real time plotting class based on Gnuplot
    '''

    def __init__(self, persist=0, debug=0):
        self.gp = Popen(['gnuplot', '-persist'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        self.plots = []

    def get_plot(self):
        self.gp = Popen(['gnuplot', '-persist'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        return self.gp

    def close_plot(self):
        self.gp.close()

    def scatter(self, x, y, names=[], title='', style='points', jitter=True):
        """
        Makes scatter plots from numpy arrays.
        if arrays are multidimensional, multiple scatter plots will be generated, pairing rows.
        :param x: array of x data points
        :param y:
        :param names:
        :param title:
        :param style:
        :param jitter:
        """
        if jitter:
            jt = numpy.random.normal(1, 1e-4, 1)[0]
        else:
            jt = 1
        if isinstance(x, numpy.ndarray):
            if not isinstance(y, numpy.ndarray):
                raise TypeError("If x is a numpy array, y must also be an array.")
            if x.shape != y.shape:
                raise ValueError("x, %s and y, %s arrays must have the same shape." % (x.shape, y.shape))
            if names:
                if len(names) != x.shape[0]:
                    raise ValueError("names list must have exactly %s items, but has %s." % (x.shape[0], len(names)))
        else:
            x = numpy.array(x)
            y = numpy.array(y)

        self.gp.stdin.write(('set title "%s"\n' % title).encode())
        if not names:
            names = ['s%s' % i for i in range(x.shape[0])]
        if len(x.shape) > 1 and len(x.shape) <= 2:
            i = 0
            for n in range(x.shape[0]):
                d = zip(x[n]*jt, y[n]*jt)
                self._plot_d(d, label=names[i], style=style)
                i += 1
            self.gp.plot(*tuple(self.plots))
        elif len(x.shape) > 2:
            pass
        else:
            d = zip(x[0]*jt, y[0]*jt)
            self._plot_d(d, label=names[0], style=style)


    def _plot_d(self, d, label="", style='points'):
        """
        Actually plots the data
        """
        self.gp.stdin.write(('set style data {}\n'.format(style)).encode())
        self.gp.stdin.flush()
        self.gp.stdin.write(("plot '-' title \"{}\" \n".format(label)).encode())
        self.gp.stdin.write(("\n".join(("%f "*len(l))%l for l in d)).encode())
        self.gp.stdin.write(b"\ne\n")
        self.gp.stdin.flush()

    def plotlines(self, data, x=None, names=[], title='', style='lines', multiplot=0):
        '''
        Create a single/multiple line plot from a numpy array or record array.
        
        :Parameters:
            - `data`: must be a numpy array or a list of lists,or record array, with series as rows
            - `x`: x values for the series: numpy array
            - `names`: is a list of strings to serve as legend labels
            - `title`: Figure Title.
            - `style`: plot styles from gnuplot: lines, boxes, points, linespoints, etc.
            - `multiplot`: Whether to make multiple subplots
        '''
        #self.gp = self.get_plot()
        if multiplot:
            sq = numpy.sqrt(len(data))
            r = numpy.floor(sq);
            c = numpy.ceil(sq)
            if len(data) == 3:
                r = 3;
                c = 1
            self.gp.stdin.write(('set multiplot layout %s,%s title "%s"\n' % (r, c, title)).encode())

        else:
            self.gp.stdin.write(('set title "%s"\n' % title).encode())
        self.gp.stdin.flush()
        if isinstance(data, list):
            data = numpy.array(data)
        if isinstance(data, numpy.core.records.recarray):
            return self._linesFromRA(data, x, style)
        if len(data.shape) > 1 and len(data.shape) <= 2:
            i = 0
            for row in data:
                if x == None:
                    x = numpy.arange(len(row))
                if names:
                    d = zip(x, row)
                    self._plot_d(d, label=names[i], style=style)
                else:
                    d = zip(x, row)
                    self._plot_d(d, style=style)
                i += 1

        elif len(data.shape) > 2:
            pass
        else:
            #            print data
            if x == None:
                x = numpy.arange(len(data))
            d = zip(x, data)
            self._plot_d(d, label=names[0], style=style)


    def _linesFromRA(self, data, x, style):
        '''
        Record-array specific code
        '''
        for n in data.dtype.names:
            if len(data.shape) > 1 and len(data.shape) <= 2:
                i = 0
                for row in data[n]:
                    if x == None:
                        x = numpy.arange(len(row))
                    d = zip(x, data)
                    self._plot_d(d, label=n + ':{}'.format(i), style=style)

                    i += 1
            elif len(data.shape) > 2:
                pass
            # TODO: figure out what to do with higher dimensional data
            else:
                d = zip(x, data[n])
                self._plot_d(d, label=n, style=style)


    def plothist(self, data, title='', names=[], multiplot=0):
        '''
        Create a sinlge/multiple Histogram plot from a numpy array or record array.
        
        :Parameters:
            - `data`: must be a numpy array or record array, with series as rows
            - `names`: is a list of strings to serve as legend labels
        '''
        #self.gp = self.get_plot()
        if multiplot:
            sq = numpy.sqrt(len(data))
            r = numpy.floor(sq);
            c = numpy.ceil(sq)
            if len(data) == 3:
                r = 3;
                c = 1
            self.stdin.write(('set multiplot layout %s,%s title "%s"\n' % (r, c, title)).encode())
        else:
            self.gp.stdin.write(('set title "%s"\n' % title).encode())
        self.gp.stdin.write(b'set style boxes\n')
        self.gp.stdin.flush()
        if isinstance(data, list):
            data = numpy.array(data)
        if isinstance(data, numpy.core.records.recarray):
            return self._histFromRA(data)
        if not names:
            names = ['s%s' % i for i in range(data.shape[0])]
        if len(data.shape) > 1 and len(data.shape) <= 2:
            for n, row in enumerate(data):
                m, bins = numpy.histogram(row, normed=True, bins=50)
                d = list(zip(bins[:-1], m))
                self._plot_d(d, label=names[n], style='boxes')


        elif len(data.shape) > 2:
            pass
        else:
            m, bins = numpy.histogram(data, normed=True, bins=50)
            d = list(zip(bins[:-1], m))
            self._plot_d(d, label=names[0], style='boxes')


    def _histFromRA(self, data):
        '''
        Record-array specific code
        '''
        for n in data.dtype.names:
            if len(data.shape) > 1 and len(data.shape) <= 2:
                i = 0
                for row in data[n]:
                    m, bins = numpy.histogram(row, normed=True, bins=50)
                    d = list(zip(bins[:-1], m))
                    self._plot_d(d, label=n + ':%s' % i, style='boxes')
                    i += 1
            elif len(data.shape) > 2:
                pass
            else:
                m, bins = numpy.histogram(data[n], normed=True, bins=50)
                d = list(zip(bins[:-1], m))
                self._plot_d(d, label=n + ':%s' % i, style='boxes')


class RTpyplot:
    """
    Creates an animated plot based on pylab
    """

    def __init__(self, nseries=1, leng=10, names=['l'], title=""):
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
            line, = P.plot(numpy.arange(self.leng), [1] * self.leng, label=self.names[i], title=self.title,
                           animated=True)
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


def start_server(server):
    server.register_instance(RTplot(persist=0))
    server.register_introspection_functions()
    server.serve_forever()


def rpc_plot(port=None):
    """
    XML RPC plot server factory function
    returns port if server successfully started or 0
    """
    if port == None:
        po = 10001
        while 1:
            if po not in __ports_used: break
            po += 1
        port = po
    try:
        server = SimpleXMLRPCServer(("localhost", port), logRequests=False, allow_none=True)
        server.register_introspection_functions()
        p = Process(target=start_server, args=(server,))
        p.daemon = True
        p.start()
    except:
        return 0
    __ports_used.append(port)
    return port

# p = Process(target=start_server)
# p.daemon = True
# p.start()
if __name__ == "__main__":
    from numpy.random import normal
    gp = RTplot()
    gp2 = RTplot()
    gp.plotlines([list(range(10)), list(range(10))], x=list(range(10, 20)), names=['a', 'b'], title='multi', style='lines', multiplot=1)
    data = normal(5,2,(3,1000))
    gp2.plothist(data=data, title="normal")
