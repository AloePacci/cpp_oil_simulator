import sys
from basic_sim.build.utils.pywrap import SIMULATOR
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
from datetime import datetime
import random
#enhanced from https://github.com/derpberk/MultiAgentEntropyDRL/blob/master/Environment/OilSpillEnvironment.py
kernel=[[1,4,7,4,1],[4,20,33,20,4],[7,33,55,33,7],[4,20,33,20,4],[1,4,7,4,1]] #for gaussian filtering
class MAP:
    def __init__(self, map_path) -> None:
        self.map=np.genfromtxt(map_path, delimiter=" ")
        self.resolution = 0
        self.visitable=[]
        for i in range(len(self.map)):
            for j in range(len(self.map[0])):
                if self.map[i][j] == 1:
                    self.visitable.append([i,j])



###Parameters
dt=0.3        #time step
kw=5    #weight of the windspeed
kc=1.5        #weight of the current speed
gamma=0.1     #weight of the random movement (brownian movement term)
flow=15     #flow of oil - Number of particles generated per time step
truerandom=True
number_of_sources=1

class SIM:
    def __init__(self,mappath=None, seed = 20) -> None:
        if mappath is None:
            self.mappath="mapas/defuniak_map.csv"
        else:
            self.mappath=mappath
        # self.mappath="mapas/ypacarai_map.csv"
        # self.mappath="mapas/algeciras.csv"
        self.map=MAP(self.mappath)
        self.map.map=self.map.map[:,:] #cpp despreciates first value
        print(f"python shape map {len(self.map.map)}, {len(self.map.map[0])}")
        self.simulator=SIMULATOR.SIMULATOR(_filepath=self.mappath, _dt=dt, _kw=kw, _kc=kc, _gamma=gamma, _flow=flow, _number_of_sources=1, _max_contamination_value=5, _source_fuel=100000, _random_seed=seed)
        self.fig = None
        self.pause=False
        self.im1=None
        self.im2=None
        self.fig_number=0
    def step(self):
        self.simulator.step()
    def reset(self):
        self.simulator.reset()

    def print_map(self):
        if self.fig_ == None:
           self.fig_=plt.figure()
        plt.clf()
        rendered = np.copy(self.map.map) * self.simulator.density
        rendered[self.map.map == 0] = np.nan
        self.im3 = plt.imshow(rendered, interpolation='nearest', cmap='jet', vmin=0, vmax=30)
        plt.scatter(self.simulator.contamination_position[1,:], self.simulator.contamination_position[0,:], c='b', label='Source points')
        plt.pause(0.0001)


    def render2(self):
        if self.im1==None:
            self.fig, self.ax = plt.subplots(1, 2, figsize=(10, 5))
            cid = self.fig.canvas.mpl_connect('key_press_event', self.onclick)
            print("press enter to Pause/resume the simulation")
            rendered = np.copy(self.map.map) * self.get_density()
            rendered[self.map.map == 0] = np.nan
            self.im1 = self.ax[1].imshow(rendered, interpolation='nearest', cmap='jet', vmin=0, vmax=30)
        else:
            self.ax[0].clear()
            self.ax[1].clear()
            rendered = np.copy(self.map.map) * self.get_density()
            rendered[self.map.map == 0] = np.nan
            self.im1 = self.ax[1].imshow(rendered, interpolation='nearest', cmap='jet', vmin=0, vmax=30)
            self.source_p = self.ax[1].scatter(self.simulator.source_points[1], self.simulator.source_points[0], c='r', label='Source points')
            plt.legend()
        # print(np.asarray(self.get_density()).max())

        self.u=self.simulator.u.tolist()
        self.v=self.simulator.v.tolist()
        for i in range(len(self.map.map)):
            for j in range(len(self.map.map[0])):
                if self.map.map[i][j]==0:
                    self.v[j][i]=0
                    self.u[j][i]=0
        self.ax[0].quiver(self.simulator.y, self.simulator.x,self.u, self.v, scale=50)
        self.ax[0].set_xlim((0, self.map.map.shape[1]))
        self.ax[0].set_ylim((0, self.map.map.shape[0]))
        self.ax[0].scatter(self.simulator.contamination_position[1], self.simulator.contamination_position[0], s=20)
        self.ax[0].invert_yaxis()

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        plt.pause(0.001)
        while self.pause and plt.fignum_exists(1):
            plt.pause(0.001)

    def onclick(self, event):
        if(event.key == "enter"):
            self.pause=not self.pause
            print("paused" if self.pause else "resumed")
        # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
        #     ('double' if event.dblclick else 'single', event.button,
        #     event.x, event.y, event.xdata, event.ydata))


    
    def savefigs_paper(self, addname=""):
        SMALL_SIZE = 8
        PICO_SIZE = 4
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 20

        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=PICO_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=PICO_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        if self.im1 ==None:
            self.im1= plt.figure(figsize=(10, 10))
            self.im1.canvas.mpl_connect('key_press_event', self.onclick)
            print("press enter to Pause/resume the simulation")
        # rendered = np.copy(self.map.map) * self.get_density()
        rendered = np.copy(self.map.map) * self.get_density(apply_gaussian=True)
        rendered_raw = np.copy(self.map.map) * self.get_density(apply_gaussian=False)
        # rendered[self.map.map == 0] = np.nan
        self.u=self.simulator.u.tolist()
        self.v=self.simulator.v.tolist()
        for i in range(len(self.map.map)):
            for j in range(len(self.map.map[0])):
                if self.map.map[i][j]==0:
                    self.v[j][i]=0
                    self.u[j][i]=0
        wind_speed=self.simulator.wind_speed.tolist()
        wx=np.copy(navigation_map.T)*wind_speed[0]
        wy=np.copy(navigation_map.T)*wind_speed[1]
        ex = np.zeros_like(navigation_map)
        
        
        cmap = cmr.get_sub_cmap('viridis', 0.60, 0.99)
        
        ################
        ##  PARTICLES
        ################
        plt.clf()
        plt.axis("off")
        plt.gca().set_xlim((0, self.map.map.shape[1]))
        plt.gca().set_ylim((0, self.map.map.shape[0]))
        plt.imshow(navigation_map*np.NaN, vmin=0, vmax=1, cmap='copper_r', aspect='auto', alpha = 1 - navigation_map, zorder=10)
        plt.imshow(ex, cmap=cmr.get_sub_cmap('viridis', 0.30, 0.99), vmin=0, vmax=1, alpha=navigation_map)
        plt.scatter(self.simulator.contamination_position[1], self.simulator.contamination_position[0], s=5, c='r', label="Particles")
        # plt.scatter(self.simulator.source_points[1], self.simulator.source_points[0], c='r', label='Source points')
        plt.imshow(navigation_map*np.NaN, vmin=0, vmax=1, cmap='copper_r', aspect='auto', alpha = 1 - navigation_map)
        plt.legend()
        plt.gca().invert_yaxis()
        plt.savefig(f'figs/particles.pdf',bbox_inches='tight', pad_inches = 0, dpi=400, transparent = True) #_{random.randint(0,2000)}

        ################
        ##  DENSITY
        ################

        plt.clf()
        plt.axis("off")
        plt.gca().set_xlim((0, self.map.map.shape[1]))
        plt.gca().set_ylim((0, self.map.map.shape[0]))
        plt.imshow(navigation_map*np.NaN, vmin=0, vmax=1, cmap='copper_r', aspect='auto', alpha = 1 - navigation_map, zorder=10)
        plt.imshow(ex, cmap=cmr.get_sub_cmap('viridis', 0.30, 0.99), vmin=0, vmax=1, alpha=navigation_map)
        plt.imshow(rendered_raw, vmin=0, vmax=1, cmap=cmap, aspect='auto', alpha=rendered)
        plt.gca().invert_yaxis()
        plt.savefig(f'figs/densityraw.pdf',bbox_inches='tight', pad_inches = 0, dpi=400, transparent = True) #_{random.randint(0,2000)}


        ################
        ##  DENSITY RAW
        ################

        plt.clf()
        plt.axis("off")
        plt.gca().set_xlim((0, self.map.map.shape[1]))
        plt.gca().set_ylim((0, self.map.map.shape[0]))
        plt.imshow(navigation_map*np.NaN, vmin=0, vmax=1, cmap='copper_r', aspect='auto', alpha = 1 - navigation_map, zorder=10)
        plt.imshow(ex, cmap=cmr.get_sub_cmap('viridis', 0.30, 0.99), vmin=0, vmax=1, alpha=navigation_map)
        plt.imshow(rendered, vmin=0, vmax=1, cmap=cmap, aspect='auto', alpha=rendered)
        plt.gca().invert_yaxis()
        plt.savefig(f'figs/density.pdf',bbox_inches='tight', pad_inches = 0, dpi=400, transparent = True) #_{random.randint(0,2000)}


        ################
        ##  CURRENTS
        ################

        plt.clf()
        plt.axis("off")
        plt.imshow(navigation_map*np.NaN, vmin=0, vmax=1, cmap='copper_r', aspect='auto', alpha = 1 - navigation_map, zorder=10)
        plt.imshow(ex, cmap=cmr.get_sub_cmap('viridis', 0.30, 0.99), vmin=0, vmax=1, alpha=navigation_map)
        # plt.imshow(rendered, interpolation='nearest', cmap='jet', vmin=0, vmax=30)
        plt.gca().set_xlim((0, self.map.map.shape[1]))
        plt.gca().set_ylim((0, self.map.map.shape[0]))
        plt.quiver(self.simulator.y[::2][::2], self.simulator.x[::2][::2],self.u[::4], self.v[::4], scale=30)
        plt.savefig(f'figs/currents.pdf') #_{random.randint(0,2000)}

        ################
        ##  WIND
        ################

        plt.clf()
        plt.axis("off")
        plt.imshow(navigation_map*np.NaN, vmin=0, vmax=1, cmap='copper_r', aspect='auto', alpha = 1 - navigation_map, zorder=10)
        plt.imshow(ex, cmap=cmr.get_sub_cmap('viridis', 0.30, 0.99), vmin=0, vmax=1, alpha=navigation_map)
        # plt.imshow(rendered, interpolation='nearest', cmap='jet', vmin=0, vmax=30)
        # plt.gca().set_xlim((0, self.map.map.shape[1]))
        # plt.gca().set_ylim((0, self.map.map.shape[0]))
        plt.quiver(self.simulator.y[::4][::2], self.simulator.x[::4][::2], wx[::8], wy[::8], scale=20)
        # plt.quiver(self.simulator.y[::8][::2], self.simulator.x[::2][::8], wx[::16], wy[::16], scale=20)
        plt.savefig(f'figs/wind.pdf') #_{random.randint(0,2000)}




    def render(self):


        if self.im1 is None:
            self.fig, self.ax = plt.subplots(1, 2, figsize=(10, 5))
            u=self.simulator.u.tolist()
            v=self.simulator.v.tolist()
            for i in range(len(self.map.map)):
                for j in range(len(self.map.map[0])):
                    if self.map.map[i][j]==0:
                        v[j][i]=0
                        u[j][i]=0
            self.ax[0].quiver(self.simulator.y, self.simulator.x, v, u, scale=50)
            self.ax[0].set_xlim((0, self.map.map.shape[1]))
            self.ax[0].set_ylim((0, self.map.map.shape[0]))
            self.ax[0].invert_yaxis()
            self.im0 = self.ax[0].scatter(self.simulator.contamination_position[1], self.simulator.contamination_position[0], s=30, lw=0.5)
            rendered = np.copy(self.map.map) * self.simulator.density
            rendered[self.map.map == 0] = np.nan
            self.im1 = self.ax[1].imshow(rendered, interpolation='nearest', cmap='jet', vmin=0, vmax=30)
            self.source_p = self.ax[1].scatter(self.simulator.source_points[1], self.simulator.source_points[0], c='r', label='Source points')
            plt.legend()

        else:
            rendered = self.simulator.density
            rendered[mask] = np.nan
            contamination_position=[self.simulator.contamination_position[1], self.simulator.contamination_position[0]]
            self.im0.set_offsets(contamination_position)
            self.im1.set_data(rendered)


        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        plt.pause(0.00001)

    def get_density(self, apply_gaussian = True):
        # if apply_gaussian:        
        #     #apply gaussian filter
        #     aux=np.zeros(self.map.map.shape)
        #     for c in self.map.visitable:
        #         i=c[0]
        #         j=c[1]
        #         found=0
        #         for a in range(5):
        #             for b in range(5):
        #                 if [i+a-2,j+b-2] in self.map.visitable:
        #                     found+=1
        #                     aux[i][j]+=self.simulator.get_normalized_density()[i+a-2][j+b-2]*kernel[a][b]
        #         aux[i][j]/=found
        #         if aux[i][j]>25:
        #             aux[i][j]=25
        #     aux[self.map.map == 0] = -1
        # else:
        aux=self.simulator.get_normalized_density(apply_gaussian)
        return aux
    

    def savefig(self, apply_gaussian = True, intensity=5, addname=""):
        if self.im1 is None:
            self.im1 = plt.figure(figsize=(10, 5))
        plt.clf()
        self.u=self.simulator.u.tolist()
        self.v=self.simulator.v.tolist()
        for i in range(len(self.map.map)):
            for j in range(len(self.map.map[0])):
                if self.map.map[i][j]==0:
                    self.v[j][i]=0
                    self.u[j][i]=0
        
        rendered = self.get_density(apply_gaussian = apply_gaussian)
    	#apply gaussian filter

        rendered[self.map.map == 0] = np.nan
        fig = plt.imshow(rendered*intensity, interpolation='nearest', cmap='gray', vmin=0, vmax=30)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        plt.scatter(self.simulator.source_points[1], self.simulator.source_points[0], c='r', label='Source points')
        plt.savefig(f'figs/{self.fig_number}grid{addname}.png', dpi=400) #_{random.randint(0,2000)}
        self.fig_number+=1


    def save_fig_evolution(self,step=0):
        if step==0:
            self.fig, self.ax = plt.subplots(1,6, figsize=(20,100/6))
            self.im1=1
        self.u=self.simulator.u.tolist()
        self.v=self.simulator.v.tolist()
        for i in range(len(self.map.map)):
            for j in range(len(self.map.map[0])):
                if self.map.map[i][j]==0:
                    self.v[j][i]=0
                    self.u[j][i]=0
        
        rendered = self.get_density(apply_gaussian = False)
    	#apply gaussian filter

        rendered[self.map.map == 0] = np.nan
        fig = self.ax[step].imshow(rendered*5, interpolation='nearest', cmap='jet', vmin=0, vmax=30)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
    def update_to_time(self, t):

        """ Update the environment to a given time """

        self.reset()

        for _ in range(t):
            self.step()

        return self.density

    def evaluate(self, pos: np.ndarray):
        """ Measure the environment at a given position """

        return self.density[pos[0].astype(int),pos[1].astype(int)]



if __name__ == "__main__":
    navigation_map = np.genfromtxt(f'mapas/defuniak_map.csv', delimiter=' ')
    start=datetime.now()
    compare=[]
    iterations = 15
    nmaps=2
    for i in range(0, nmaps):
        # np.random.seed(i)
        a=SIM("mapas/defuniak_map.csv", seed = 10)
        print(f"new map {i}")
        a.reset()
        step=0
        for _ in range(iterations+1):
            # if _%(iterations//5) == 0:
            #     a.save_fig_evolution(step)
            #     step+=1
            #     print(f"particle number: {len(a.simulator.contamination_position[0])}")
            a.step()
            if _ == 50:
                a.savefigs_paper()
    plt.savefig(f'figs/grid.png', dpi=400) #_{random.randint(0,2000)}
    afterquery=datetime.now()
    difference = afterquery - start
    
    aux=divmod(difference.seconds, 60)
    print(f"query of {iterations} delayed {difference.seconds}.{difference.microseconds} seconds")
    # plt.show()
