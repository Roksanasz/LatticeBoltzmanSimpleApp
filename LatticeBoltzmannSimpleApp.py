from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import sys, time, random
from matplotlib.pyplot import imshow

import numpy as np


class WorkerObject(QtCore.QObject):

    signalStatus = QtCore.pyqtSignal(object)

    def __init__(self, parent=None):
        super(self.__class__, self).__init__(parent)
        #parameters
        Nx                     = 1600    # resolution x-dir
        Ny                     = 600    # resolution y-dir
        rho0                   = 100    # average density
        
        self.Nt                = 400   # number of timesteps
        self.plotRealTime = True # switch on for plotting as the simulation goes along
        # Lattice speeds / weights
        NL = 9
        idxs = np.arange(NL)
        cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
        cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
        inv = np.array([0,3,4,1,2,7,8,5,6])
        weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) # sums to 1
        
        # Initial Conditions
        self.F = np.ones((Ny,Nx,NL))
        np.random.seed(42)
        self.F += 0.01*np.random.randn(Ny,Nx,NL)
        X, Y = np.meshgrid(range(Nx), range(Ny))
        self.F[:,:,3] += 1 * (1+0.2*np.cos(2*np.pi*X/Nx*4))
        rho = np.sum(self.F,2)
        for i in idxs:
            self.F[Ny-1,Nx-1,i] = 0
        
        for i, cx, cy in zip(idxs, cxs, cys):
            self.F[:,:,i] = np.roll(self.F[:,:,i], cx, axis=1)
            self.F[:,:,i] = np.roll(self.F[:,:,i], cy, axis=0)
              
        Nx = 1600
        Ny = 600
        self.data = np.zeros((Ny,Nx))
        self.x = np.arange(0,Nx)
        self.y = np.arange(0,Ny)
        self.Y, self.X = np.meshgrid(self.x, self.y)
        #ux = 0.6
        #uy = 0.1
        cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
        cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
        self.f = np.zeros(self.F.shape)
        for i, cx, cy, w in zip(idxs, cxs, cys, weights):
            self.f = lambda ux, uy: rho0 * w * ( 1 - 3/2*(ux*ux+uy*uy)  + 3*(cx*ux+cy*uy) + 9/2*(cx*ux+cy*uy)*(cx*ux+cy*uy))
        

    @QtCore.pyqtSlot()        
    def startWork(self):
        print("StartWork")
        x0 = 250
        y0 = 400
        
        while 1 > 0:
            x0, y0 = self.doWork(x0, y0)
            

    def doWork(self, x0, y0):         
        NL = 9    
        idxs = np.arange(NL)
        cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
        cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
        Nx = 1600
        Ny = 600
        self.F = np.ones((Ny,Nx,NL))
        np.random.seed(42)
        self.F += 0.01*np.random.randn(Ny,Nx,NL)
        self.Y, self.X = np.meshgrid(self.x, self.y)
        self.F[:,:,3] += 1 * (1+0.2*np.cos(2*np.pi*self.X/Nx*4))
        rho = np.sum(self.F,2)
        for i in idxs:
            self.F[Ny-1,Nx-1,i] = 0
        
       
        # Drift
        for i, cx, cy in zip(idxs, cxs, cys):
            self.F[:,:,i] = np.roll(self.F[:,:,i], cx, axis=1)
            self.F[:,:,i] = np.roll(self.F[:,:,i], cy, axis=0)
# przeszkody
        #od lewej trzecia przeszkoda w kształcie koła
        X, Y = np.meshgrid(range(Nx), range(Ny))
        shape3 = (X - Nx/2)**2 + (Y - Ny/2)**2 < (Ny/10)**2
        #od lewej pierwsza przeszkoda
        X, Y = np.meshgrid(range(Nx), range(Ny))
        shape = (X - Nx/22)**2 + (Y - Ny/2)**2 < (Ny/10)**2
        #od lewej środkowa przeszkoda
        X, Y = np.meshgrid(range(Nx), range(Ny))
        shape2 = (X - Nx/3)**2 + (Y - Ny/8)**2 < (Ny/10)**2
        #ścianka dolna
        X, Y = np.meshgrid(range(Nx), range(Ny))
        shape4 = Y  < (Ny/70)
        #ścianka górna
        X, Y = np.meshgrid(range(Nx), range(Ny))
        shape5 = (Y - Ny/1)**2 < (Ny/60)**2
        inv = np.array([0,3,4,1,2,7,8,5,6])
        bndryF = self.F[shape,:]
        bndryF = bndryF[:,inv]
        bndryF2 = self.F[shape2,:]
        bndryF2 = bndryF2[:,inv]
        bndryF3 = self.F[shape3,:]
        bndryF3 = bndryF3[:,inv]
        bndryF4 = self.F[shape4,:]
        bndryF4 = bndryF4[:,inv]
        bndryF5 = self.F[shape5,:]
        bndryF5 = bndryF5[:,inv]        
        NL = 9
        tau                    = 1   # collision timescale
        weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) # sums to 1
        #cx = 0.6
        #cy = 0.2
        #dx, dy = random.randint(0,1)*2-1, random.randint(0,1)*2-1
        rho = np.sum(self.F,2)
        ux  = np.sum(self.F*cxs,2) / rho
        uy  = np.sum(self.F*cys,2) / rho
        for i, cx, cy, w in zip(idxs, ux, uy, weights):
            self.data[:,:] = (self.f(ux,uy)-self.F[:,:,i])
        #-self.F[:,:,i]
        self.F[:,:,i] += -(1.0/tau) * (self.f(ux,uy)-self.F[:,:,i])
        ## Apply boundary 
        self.F[shape,:] = bndryF
        self.F[shape2,:] = bndryF2
        self.F[shape3,:] = bndryF3
        self.F[shape4,:] = bndryF4
        self.F[shape5,:] = bndryF5
        #Nt                     = 400
        #if (True):
            
        ux[shape] = 0
        uy[shape] = 0
        ux[shape2] = 0
        uy[shape2] = 0
        ux[shape3] = 0
        uy[shape3] = 0
        ux[shape4] = 0
        uy[shape4] = 0
        ux[shape5] = 0
        uy[shape5] = 0

        self.vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) - (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
        self.vorticity[shape] = np.nan
        self.vorticity[shape2] = np.nan
        self.vorticity[shape3] = np.nan
        self.vorticity[shape4] = np.nan
        self.vorticity[shape5] = np.nan
        #cmap.set_bad('black')
        self.signalStatus.emit(self.vorticity)
        
        time.sleep(0.1)
        return  cx, cy

class App(QtWidgets.QMainWindow):

    signalStatus = QtCore.pyqtSignal(object)
    abortSignal = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super(App, self).__init__(parent)

        self.button_start = QtWidgets.QPushButton('Start',self)
        self.button_start.setFixedWidth(120)
        self.button_start.setStyleSheet(#setting variable margins
        "*{margin-left: " + str(10) +"px;"+
        "margin-right: " + str(10) +"px;"+
        '''
        border: 0.5em solid '#264653';
        color: Black;
        font-family: 'shanti';
        font-size: 6em;
        border-radius: 2em;
        padding: 2em 0;
        margin-top: 0.3em;
        }
        *:hover{
            background: '#264653';
            
        }
        ''')
        self.button_cancel = QtWidgets.QPushButton('Stop', self)
        self.button_cancel.setFixedWidth(120)
        self.button_cancel.setStyleSheet(#setting variable margins
        "*{margin-left: " + str(10) +"px;"+
        "margin-right: " + str(10) +"px;"+
        '''
        border: 0.5em solid '#E76F51';
        color: Black;
        font-family: 'shanti';
        font-size: 6em;
        border-radius: 2em;
        padding: 2em 0;
        margin-top: 0.3em;
        }
        *:hover{
            background: '#E76F51';
            
        }
        ''')
        self.label_status = QtWidgets.QLabel('', self)

        self.mainbox = QtWidgets.QWidget(self)
        self.layout = QtWidgets.QVBoxLayout()
        self.setWindowTitle("LBM1")
        self.mainbox.setLayout(self.layout)
        self.mainbox.setStyleSheet("QWidget { background-color: %s }" % QtGui.QColor(233, 196, 106, 255).name())
        self.setCentralWidget(self.mainbox)
        self.layout.addWidget(self.button_start)
        self.layout.addWidget(self.button_cancel)
        self.layout.addWidget(self.label_status)

        self.fig = Figure((4.0,2.0), dpi=80)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.canvastoolbar = NavigationToolbar(self.canvas, self)
        self.fig.patch.set_alpha(0.0)
        self.ax = self.fig.add_subplot(111)
        self.x =  np.arange( 0,1600 ); self.y = np.arange( 0,600 )
        self.im = self.ax.imshow(np.zeros((600, 1600)),  cmap='jet', vmin=0, vmax=1 )
        #self.pv,  = self.ax.plot( np.zeros(600) ,self.y , color="white" , alpha=0.6, lw=2 )
        #self.ph,  = self.ax.plot( self.x ,np.zeros(1600) , color="white" , alpha=0.6, lw=2)
        self.ax.set_xlim([0,1600]); self.ax.set_ylim([0,600])
        self.ax.get_yaxis().set_visible(False)
        self.ax.get_xaxis().set_visible(False)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.canvastoolbar)

        self.initWorker()



    def initWorker(self):
        self.worker = WorkerObject()
        self.worker_thread = QtCore.QThread()
        self._connectSignals()
        self.worker.moveToThread(self.worker_thread)
        self.worker_thread.start()


    def _connectSignals(self):
        self.button_start.clicked.connect(self.worker.startWork)
        self.button_cancel.clicked.connect(self.forceWorkerQuit)
        self.worker.signalStatus.connect(self.updateStatus)


    def forceWorkerQuit(self):
        print("calculation aborted")
        if self.worker_thread.isRunning():
            self.worker_thread.terminate()
        self.worker_thread.start()


    @QtCore.pyqtSlot(object)
    def updateStatus(self, obj):
        self.im.set_data(obj)
        argm = np.unravel_index(np.argmax(obj), (600,1600))
        #self.pv.set_data(obj[:,argm[0]]*250, self.y)
        #self.ph.set_data(self.x, obj[argm[1], :]*250)
        self.fig.canvas.draw()


if __name__=='__main__':
    app = QtWidgets.QApplication(sys.argv)
    thisapp = App()
    thisapp.show()
    sys.exit(app.exec_())