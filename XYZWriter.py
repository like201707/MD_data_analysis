import Lammpstrj_reader
import numpy as np
import sys
import time


def count_atom(A_type):
    # count the number of A atom in each frame
    global nve
    Atoms = 0
    for j in range(nve.nAtoms):
        if nve.atomT[0][j][0] == A_type:
            Atoms += 1
    return Atoms


def get_coordinates_velovoties(A_type):
    # get specific type atom coordinates from the trjectory
    global nve
    Atoms = count_atom(A_type)
    A_atomC = np.empty([nve.nFrames, Atoms, 3], dtype=np.float)
    A_atomV = np.empty([nve.nFrames, Atoms, 3], dtype=np.float)
    for frame in range(nve.nFrames):
        i = 0
        for atom in range(nve.nAtoms):
            if nve.atomT[frame][atom][0] == A_type:
                A_atomC[frame][i] = nve.atomC[frame][atom]
                A_atomV[frame][i] = nve.atomV[frame][atom]
                i += 1
    return (A_atomC, A_atomV)


def write_file(A_type):
    global nve
    print ("\n\nWriting xyz file ...")
    starttime = time.time()
    A_atomC = get_coordinates_velovoties(A_type)[0]
    A_atomV = get_coordinates_velovoties(A_type)[1]
    A_num = count_atom(A_type)
    with open("Type"+str(A_type)+".xyz", "w") as g:
        for i in range(nve.nFrames):
            g.write(str(A_num)+"\n")
            g.write("Frame number:"+" "+str(i+1)+"\n")
            for j in range(A_num):
                g.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(
                        A_type, A_atomC[i][j][0], A_atomC[i][j][1],
                        A_atomC[i][j][2], A_atomV[i][j][0],
                        A_atomV[i][j][1], A_atomV[i][j][2]))
        print ("Done! ..............\n\n")
        endtime = time.time()
        print ('Writing xyz file total time:%f') % (endtime-starttime)


def read_file():
    global nve
    args = sys.argv
    filename = args[1]
    nve = Lammpstrj_reader.Lammpstrj_reader(filename)


if __name__ == "__main__":
    read_file()
    write_file(13)
    write_file(9)
