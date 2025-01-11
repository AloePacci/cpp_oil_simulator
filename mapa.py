import numpy as np

class MAP:
    def __init__(self, map_path) -> None:
        self.map=np.genfromtxt(map_path, delimiter=" ")
        self.resolution = 0
        self.visitable=[]
        for i in range(len(self.map)):
            for j in range(len(self.map[0])):
                if self.map[i][j] == 1:
                    self.visitable.append([i,j])


