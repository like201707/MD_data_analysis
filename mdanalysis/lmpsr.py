# -*- coding: utf-8 -*-

import numpy as np
from progress.bar import FillingCirclesBar 

class LMPS_Reader(object):
    # format lammps dump trjectory
    def __init__(self, filename):
        
        self.filename = filename
        self.nAtoms = self.n_atoms()
        self.nFrames = self.n_frame()
        self.file = open(self.filename, "r")
        self.atomT = np.empty([self.nFrames, self.nAtoms], dtype=np.int)
        self.atomC = np.empty([self.nFrames, self.nAtoms, 3], dtype=np.float)
        self.atomV = np.empty([self.nFrames, self.nAtoms, 3], dtype=np.float)
        self.boxSize = np.empty([self.nFrames, 3], dtype=np.float)
        self.read_all_frames()

    def n_atoms(self):
        with open(self.filename) as f:
            for i in range(3):
                f.readline()
            n = f.readline()
            return int(n)

    def n_frame(self):
        linesPerFrame = self.nAtoms + 9
        count = 0
        offsets = list()
        with open(self.filename) as f:
            line = True
            while line:
                if not count % linesPerFrame:
                    offsets.append(f.tell())
                line = f.readline()
                count += 1
        n_frames = int(count/linesPerFrame)
        self.offsets = offsets
        return n_frames
        

    def read_all_frames(self):
        bar = FillingCirclesBar('Processing', max=self.nFrames)
        for frame in range(self.nFrames):
            self.file.seek(self.offsets[frame])
            self.read_one_frame(frame)
            bar.next()
        bar.finish()

    def read_one_frame(self, frame):
        f = self.file
        for i in range(5):
            f.readline()
        for i in range(3):
            line = f.readline().split()
            self.boxSize[frame][i] = float(line[1]) - float(line[0])
        f.readline()
        for i in range(self.nAtoms):
            lineinfo = f.readline().split()
            self.atomT[frame][i] = int(lineinfo[2])
            self.atomC[frame][i] = np.array(map(float, lineinfo[4:7]),
                                            dtype=np.float)
            self.atomV[frame][i] = np.array(map(float, lineinfo[7:10]),
                                            dtype=np.float)

    def close(self):
        """Close trajectory file if it was open."""
        if self.file is None:
            return
        self.file.close()
        
        
        
        
        
        
        
        
