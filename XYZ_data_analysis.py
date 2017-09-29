import XYZReader
import numpy as np
import pylab as pl
import sys
import math
import time


class Diffusion_coefficent(object):

    def __init__(self, atomC, atomV):
        self.c = atomC
        self.v = atomV
        self.nFrames, self.nAtoms, random = self.v.shape
        print ("------------------------------------")
        print ("--calculating diffusion coefficent--")
        print ("------------------------------------")
        print ('\n')
        self.max_t = 1000

    def v_auto_correlation(self):
        # velocity veloccity autocorrelation function
        # one way to calculate the diffusion coefficent
        print ("--velocity veloccity autocorrelation function--")
        self.ctt = np.empty([self.max_t], dtype=float)
        for t in range(self.max_t):
            print "\r", "--calulating progress--" + " " + \
                   str(t/10) + "%",
            ct = 0.0
            for i in range(self.nFrames):
                if i + t < self.nFrames:
                    for j in range(self.nAtoms):
                        ct += np.dot(self.v[i][j], self.v[i+t][j])
            ct = ct / (self.nAtoms * (self.nFrames - t))
            self.ctt[t] = ct
        return self.ctt

    def vv_diffusion(self):
        d = 1.0 / 3.0 * np.trapz(self.ctt)
        return d

    def root_mean_square_displacement(self):
        # mean square displacement
        # another way to calculate the diffusionc coefficent
        print ("--mean square displacement--")
        self.msd = np.empty([self.max_t], dtype=float)
        for t in range(self.max_t):
            print "\r", "--calulating progress--" + " " + \
                   str(t/10) + "%",
            msd = 0.0
            mid = 0.0
            for i in range(self.nFrames):
                if i + t < self.nFrames:
                    for j in range(self.nAtoms):
                        mid += np.square(np.subtract(self.c[i][j],
                                         self.c[i+t][j]))
                        msd = np.sum(mid)
            self.msd[t] = msd / (self.nAtoms * (self.nFrames - t))
        return self.msd

    def plot_vv(self):
        X = np.arange(0, 1000, 1)
        Y = self.ctt
        pl.plot(X, Y, 'g-', linewidth=2)
        pl.xlabel("t(fs)", size=20)
        pl.ylabel(r"$<V(0)\cdot V(t)>$", size=20)
        pl.xticks(size=20)
        pl.yticks(size=20)
        pl.show()

    def plot_msd(self):
        X = np.arange(0, 1000, 1)
        Y = self.msd
        pl.plot(X, Y, 'g-', linewidth=2)
        pl.xlabel("t(fs)", size=20)
        pl.ylabel("MSD")
        pl.xticks(size=20)
        pl.yticks(size=20)
        pl.show()


class Hopping(object):
    # The secquence of passing these two parameters
    # determin who is the central atom
    def __init__(self, A_atomC, B_atomC):
        self.A_atomC = A_atomC
        self.B_atomC = B_atomC
        self.nFrames, self.nAtoms, random = self.A_atomC.shape
        self.x_size = 41.22
        self.x_cutoff = 4.35

    def distance_A_B(self):
        # Get the distance between atom A and atom B
        # the secquence determin which atom is the center atom
        # In this case atom A is the central atom
        d_A_B = np.empty([self.nFrames, self.nAtoms, self.nAtoms],
                         dtype=np.float)
        for i in range(self.nFrames):
            for j in range(self.nAtoms):
                for k in range(self.nAtoms):
                    dx = self.A_atomC[i][j][0] - self.B_atomC[i][k][0]
                    dy = self.A_atomC[i][j][1] - self.B_atomC[i][k][1]
                    dz = self.A_atomC[i][j][2] - self.B_atomC[i][k][2]
                    if dx > self.x_size * 0.5:
                        dx -= self.x_size
                    elif dx <= -self.x_size * 0.5:
                        dx += self.x_size
                    if dy > self.x_size * 0.5:
                        dy -= self.x_size
                    elif dy <= -self.x_size * 0.5:
                        dy += self.x_size
                    if dz > self.x_size * 0.5:
                        dz -= self.x_size
                    elif dz <= -self.x_size * 0.5:
                        dz += self.x_size
                    d_A_B[i][j][k] = (dx**2 + dy**2 + dz**2)**0.5
        return d_A_B

    def help_array(self):
        d_A_B = self.distance_A_B()
        a = np.array([[[0 for k in xrange(self.nAtoms)] for j
                     in xrange(self.nAtoms)] for i in xrange(self.nFrames)])
        for i in range(self.nFrames):
            for j in range(self.nAtoms):
                for k in range(self.nAtoms):
                    if d_A_B[i][j][k] < self.x_cutoff:
                        a[i][j][k] = 1
                    else:
                        a[i][j][k] = 0
        return a

    def hopping_array_specific(self, frame_start, frame_end):
        # for each S atom
        a = self.help_array()
        hopping_times = list()
        for j in range(self.nAtoms):
            H = 0
        # for each Na atom
            for k in range(self.nAtoms):
                Flag = False
                for i in range(frame_start, frame_end):
                    if a[i][j][k] == 1 and Flag:
                        H += -1
                        Flag = False
                    if a[i][j][k] != a[i+1][j][k]:
                        if a[i][j][k] == 1:
                            H += 1
                            Flag = True
            hopping_times.append(H)
        return hopping_times

    def hopping_rate(self, frame_start, frame_end):
        counter = 0
        hopping_times = self.hopping_array(frame_start, frame_end)
        for num in hopping_times:
            total_times = sum(hopping_times)
            if num != 0:
                counter += 1
        return (counter/float(self.nAtoms), total_times)

    def flag(self):
        b = np.empty([self.nAtoms, self.nAtoms], dtype=bool)
        for j in range(self.nAtoms):
            for k in range(self.nAtoms):
                b[j][k] = False
        return b

    def hopping_time(self, frame_start, frame_end):
        # for each S atom
        b = self.flag()
        a = self.help_array()
        hopping_times = list()
        H = 0
        for i in range(frame_start, frame_end):
            for j in range(self.nAtoms):
                for k in range(self.nAtoms):
                    if a[i][j][k] == 1 and b[j][k]:
                        H += -1
                        b[j][k] = False

                    if a[i][j][k] != a[i+1][j][k]:
                        if a[i][j][k] == 1:
                            H += 1
                            b[j][k] = True
            hopping_times.append(H)
        return hopping_times



class Ionpair_probability(Hopping):
    def counter(self):
        # A is the central atom
        C_A_B = np.empty([self.nFrames, self.nAtoms], dtype=np.int)
        d_A_B = self.distance_A_B()
        for i in range(self.nFrames):
            for j in range(self.nAtoms):
                counter_A_B = 0
                for k in range(self.nAtoms):
                    if d_A_B[i][j][k] < self.x_cutoff:
                        counter_A_B += 1
                C_A_B[i][j] = counter_A_B
        return C_A_B

    def probability(self):
        C_A_B = self.counter()
        prb = dict()
        for key in range(5):
            counter = 0
            for i in range(self.nFrames):
                for j in range(self.nAtoms):
                    if C_A_B[i][j] == key:
                        counter += 1
            value = counter/float(self.nAtoms * self.nFrames)
            prb[key] = value
        return prb


class R_dis_funtion(Hopping):
    def gofr(self):
        self.dr = 0.1
        self.pi = math.pi
        self.bulk_d = self.nAtoms/(self.x_size ** 3)
        self.rmax = 10
        self.n = int(self.rmax/self.dr)
        d_A_B = self.distance_A_B()
        r = np.zeros([self.n], dtype=np.int)
        g = np.empty([self.n], dtype=np.float)
        c = np.zeros([self.n], dtype=np.float)
        for i in range(self.nFrames):
            for j in range(self.nAtoms):
                for k in range(self.nAtoms):
                    if int(d_A_B[i][j][k]) < self.rmax:
                        r[int(d_A_B[i][j][k]/self.dr)] += 1
        constant = self.nAtoms * self.nFrames * self.bulk_d
        c_help = 0
        for n in range(self.n):
            g[n] = r[n]/(constant*4 * self.pi * self.dr * (self.dr*(n+1))**2)
            c_help += r[n]/float((self.nAtoms * self.nFrames))
            c[n] = c_help
        return (g, c)

    def plot(self):
        tmp = self.gofr()
        g = tmp[0]
        c = tmp[1]
        x = np.array([0.05 + self.dr * i for i in range(self.n)])
        fig, ax1 = pl.subplots()
        ax2 = ax1.twinx()
        ax1.plot(x, g, 'g-', linewidth=2)
        ax1.set_xlabel(r'$r(Na-S)(\AA)$', fontsize=20)
        ax1.set_ylabel('g(r)', color='k', fontsize=20)
        ax1.xaxis.set_tick_params(labelsize=20)
        ax1.yaxis.set_tick_params(labelsize=20)
        ax2.plot(x, c, 'g--', linewidth=2)
        ax2.set_ylabel('n(r)', color='k', fontsize=20)
        ax2.yaxis.set_tick_params(labelsize=20)
        pl.xlim([0, 6])
        fig.tight_layout()
        pl.show()


def read_file():
    global A, B
    args = sys.argv
    filenameA = args[1]
    A = XYZReader.XYZReader(filenameA)
    try:
        filenameB = args[2]
        B = XYZReader.XYZReader(filenameB)
    except IndexError:
        pass


def hopping_cal():
    global A, B
    read_file()
    TASK1 = Hopping(A.atomC, B.atomC)
    print (TASK1.hopping_array_specific(0, 9999))
    print (TASK1.hopping_time(0,9999))
#    print (TASK1.hopping_array(0, 9999), TASK1.hopping_rate(0, 9999))
#    print (TASK1.hopping_array(10000, 19999), TASK1.hopping_rate(10000, 19999))
#    print (TASK1.hopping_array(20000, 29999), TASK1.hopping_rate(20000, 29999))
#    print (TASK1.hopping_array(0, 29999), TASK1.hopping_rate(0, 29999))
#    TASK1.hopping_array_avg(7000, 8000)

def ionpair_cal():
    global A, B
    read_file()
    TASK2 = Ionpair_probability(A.atomC, B.atomC)
    print (TASK2.probability())


def gofr_cal():
    global A, B
    start = time.time()
    read_file()
    TASK3 = R_dis_funtion(A.atomC, B.atomC)
    print (TASK3.gofr()[0])
    print ("-----------------------------")
    print (TASK3.gofr()[1])
    TASK3.plot()
    end = time.time()
    print (end - start)


def vv_cal():
    global A
    read_file()
    start = time.time()
    TASK4_A = Diffusion_coefficent(A.atomC, A.atomV)
    print (TASK4_A.v_auto_correlation())
    print (TASK4_A.vv_diffusion())
    TASK4_A.plot_vv()
    end = time.time()
    print ("Total time: {}".format(end - start))


def rmsd_cal():
    global A
    read_file()
    start = time.time()
    TASK4_B = Diffusion_coefficent(A.atomC, A.atomV)
    print (TASK4_B.root_mean_square_displacement())
    TASK4_B.plot_msd()
    end = time.time()
    print ("Total time: {}".format(end - start))


if __name__ == "__main__":
    hopping_cal()
