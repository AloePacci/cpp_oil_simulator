import torch
from torch.utils.data.sampler import BatchSampler, SubsetRandomSampler
from torch.utils.data import Dataset
from torchvision import datasets
from torchvision.transforms import ToTensor
from basic_sim_fusion import SIM
import numpy as np
import sys


class MapDataset(Dataset):
    def __init__(self, map_simulation_depth, num_different_maps, memory_steps=10, map_interstep=10, map_name="acoruna_port" , map_path="cpp_py_oil_simulator/mapas/", load=False, transform=None):
        self.mapa=f"{map_path}{map_name}.csv"
        self.transform=transform
        self.contamination_positions_=[]
        map_simulation_depth += 1 #add one more as last data needs next for prediction 
        
        print(f"Creating Dataset for map {self.mapa}")
        if load: #TODO: implement load
            return
        else:
            #charge simulator
            a=SIM(self.mapa)
            for i in range(num_different_maps):
                print(f"simulating map {i+1} of {num_different_maps}. ", end="") #TODO: create a situation where there is no contamination
                a.reset()
                start_offset=sum(sum(a.get_density(apply_gaussian=False)))
                print(f"map started with contamination at point [{a.simulator.source_points[1]},{a.simulator.source_points[0]}]")
                #initialize map and data array to memory_steps
                data=[]
                for i in range(memory_steps*map_interstep):
                    a.step()
                    if i%map_interstep == 0:
                        data.append(np.array(a.get_density(apply_gaussian=False)))
                    
                self.contamination_positions_.append(np.array(data))
                for _ in range(0,map_simulation_depth*map_interstep): #now append after each new data
                    if _%map_interstep == 0:
                        data.pop(0)
                        data.append(np.array(a.get_density(apply_gaussian=False)))
                        # print((sum(sum(data))))
                        self.contamination_positions_.append(np.array(data))
                        if _*100%(map_simulation_depth*map_interstep)==0:
                            print(f"process at {100*_/(map_simulation_depth*map_interstep)}%", end='\r')
                    a.step()
                print(f"map simulation ended with {len(a.simulator.contamination_position[0])} particles, density {sum(sum(self.contamination_positions_[-1][-1]))-start_offset}, and {map_simulation_depth} episodes")
            # print(np.sum(np.sum(self.contamination_positions_)))
            contamination_positions=np.array(self.contamination_positions_)
            # print(sum(sum(contamination_positions))-start_offset)
            self.contamination_positions=torch.from_numpy(contamination_positions)
            # print(sum(sum(self.contamination_positions)-start_offset))
            # print(f"new simulation has {sum(sum(sum([x for x in self.contamination_positions.numpy()])))}")
            
    def __len__(self):
        return len(self.contamination_positions)-1

    def to(self, device):
        self.contamination_positions = self.contamination_positions.to(device)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
            
        landmarks = self.contamination_positions[idx]
        future = self.contamination_positions[idx+1][-1]
        sample = {'density': landmarks, 'future': future}

        if self.transform:
            sample = self.transform(sample)

        return sample

if __name__ == "__main__":
    a=MapDataset(map_simulation_depth=5, num_different_maps=10, memory_steps=10, map_interstep=10, map_name="acoruna_port" , map_path="mapas/", load=False, transform=None)