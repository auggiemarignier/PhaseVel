from PxMCMC import Params, PxMCMC
import numpy as np
import matplotlib.pyplot as plt

# Simple checkerboard example

class Checkerboard:
    def __init__(self,rows=8,cols=8,blocksize=3):
        self.rows = rows
        self.cols = cols
        self.blocksize = blocksize

        self.fill_board()

    def fill_board(self):
        self.board = np.kron([[1, 0] * (self.cols//2), [0, 1] * (self.cols//2)] * (self.rows//2), np.ones((self.blocksize, self.blocksize)))


    def plot_board(self):
        plt.imshow(self.board)
        plt.show()
                


class Simple(PxMCMC):
    '''
    Need to override forward, and gradg to ignore ylmavs
    '''
    pass


if __name__ == '__main__':
    board = Checkerboard()
    board.plot_board()