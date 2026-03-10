

import numpy as np
import pyqtgraph as pg
from PyQt6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from PyQt6.QtCore import QTimer
from typing import Callable

FloatGetter = Callable[[], float]
FloatSetter = Callable[[float], None]

class Oscilloscope(QMainWindow):
    def __init__(self, getter: FloatGetter, buffer_size: int = 1000):
        super().__init__()

        self.y_dense = np.zeros(buffer_size)

        self.getter = getter

        self.dt = 0.01  #to powinno byc polaczone z dt  ATAT

        self.win = pg.GraphicsLayoutWidget(title="title1")
        self.setCentralWidget(self.win)

        self.plot = self.win.addPlot(title="title2")
        self.plot.setYRange(-5.0, 5.0)
        self.plot.addLegend()
        self.plot.showGrid(x=True, y=True)

        self.curve_y = self.plot.plot(pen='y', name="output")
        #self.curve_u = self.plot.plot(pen=pg.mkPen('g', style=pg.QtCore.Qt.PenStyl.CashLine), name="input")

        self.timer = QTimer()
        self.timer.timeout.connect(self.update)

    def update(self):
        y_new = self.getter()
        self.y_dense = np.roll(self.y_dense, -1)     
        self.y_dense[-1] = y_new
        self.curve_y.setData(self.y_dense)
        self.timer.start(16)


'''
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Oscilloscope()
    window.resize(800, 400)
    window.show()
    sys.exit(app.exec()) # tu jest petla while
'''