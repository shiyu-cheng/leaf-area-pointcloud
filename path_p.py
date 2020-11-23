import math
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import numpy as np
from numpy import arange
import bisect


class Windows(object):
    """
    divide the space into windows and calculate the stat of light beams in windows
    unit: degree
    """
    def __init__(self):

        self.resolution = 0.1
        self.n_bins = 25

        self.width = 0.5
        self.height = 0.5
        # assume the step is equal or larger than the window side length
        self.width_step = self.width
        self.height_step = self.height
        self.is_gap_threshold = 5

        self.windows = {}
        self.windows_is_gap = {}
        self.start_zeniths = []
        self.start_azimuths = []

        for az in arange(0, 360, self.width_step):
            self.start_azimuths.append(az)
        
        for z in arange(0, 90, self.height_step):
            self.start_zeniths.append(z)
            self.windows[z] = {}
            self.windows_is_gap[z] = {}
            for az in arange(0, 360, self.width_step):
                self.windows[z][az] = []
                self.windows_is_gap[z][az] = False

    def add_to_window(self, z, az):
        # if int(z * 10300501) % 10000 == 0:
        #     print(str(z) + ' ' + str(az))
        index = bisect.bisect(self.start_zeniths, z)
        # start_z = the first zenith smaller than z
        start_z = self.start_zeniths[index-1]
        # if z is in the window
        if z < start_z + self.height:
            index = bisect.bisect(self.start_azimuths, az)
            start_az = self.start_azimuths[index-1]
            if az < start_az + self.width:
                self.windows[start_z][start_az].append(az)

    def get_slice_boundary(self, slice_start_z, slice_step):
        index = bisect.bisect(self.start_zeniths, slice_start_z)
        start_z = self.start_zeniths[index-1]
        index = bisect.bisect(self.start_zeniths, slice_start_z + slice_step)
        end_z = self.start_zeniths[index-1]
        return start_z, end_z
    
    def get_gap_histogram(self, start_z, end_z):
        ray_cnt = self.width * self.height / self.resolution / self.resolution
        histogram = np.zeros(self.n_bins)
        all_gap_cnt = 0
        all_ray_cnt = 0
        max_neg_log_gap_p = 0
        sample_cnt = 0
        big_gap_window_cnt = 0
        all_window_cnt = 0
        # get l_max (which is the max -log(p) value)
        for z in arange(start_z, end_z, self.height_step):
            for az in arange(0, 360, self.width_step):
                all_window_cnt += 1
                if self.windows_is_gap[z][az] == False:
                    sample_cnt += 1
                    non_gap_cnt = len(self.windows[z][az])
                    # if there is no gap in the window, add half ray as a gap
                    if non_gap_cnt == ray_cnt:
                        non_gap_cnt = ray_cnt - 0.5
                    neg_log_gap_p = -math.log(1 - non_gap_cnt / ray_cnt)
                    if max_neg_log_gap_p < neg_log_gap_p:
                        max_neg_log_gap_p = neg_log_gap_p * 1.0001
                else:
                    big_gap_window_cnt += 1
        # populate the histogram
        for z in arange(start_z, end_z, self.height_step):
            for az in arange(0, 360, self.width_step):
                if self.windows_is_gap[z][az] == False:
                    non_gap_cnt = len(self.windows[z][az])
                    if non_gap_cnt == ray_cnt:
                        non_gap_cnt = ray_cnt - 0.5
                    neg_log_gap_p = -math.log(1 - non_gap_cnt / ray_cnt)
                    histogram[int(neg_log_gap_p/max_neg_log_gap_p * self.n_bins)] += 1 / sample_cnt
                    all_gap_cnt += ray_cnt - non_gap_cnt
                    all_ray_cnt += ray_cnt
        '''for i in arange(len(histogram)):
            print(histogram[i], '')
        print('')'''
        if all_ray_cnt != 0:
            overall_gap_p = all_gap_cnt / all_ray_cnt
        else:
            overall_gap_p = 1
        return overall_gap_p, histogram, big_gap_window_cnt / all_window_cnt

    def tag_gap(self):
        for z in arange(0, 90, self.height_step):
            az_base = -1
            az_static = 0
            for az in arange(0, 360, self.width_step):
                az_static = az
                # register the start of a list of gap
                if len(self.windows[z][az]) == 0:
                    if az_base == -1:
                        az_base = az
                # the end of a lists of gap
                elif az_base != -1: 
                    if az - az_base >= self.is_gap_threshold:
                        for i in arange(az_base, az, self.width_step):
                            self.windows_is_gap[z][i] = True
                    az_base = -1
            if az_base != -1:
                for az in arange(360, 720, self.width_step):
                    if len(self.windows[z][az-360]) != 0:
                        break
                    az_static = az
                if az_static - az_base >= self.is_gap_threshold:
                    for az in arange(az_base, 360, self.width_step):
                        self.windows_is_gap[z][az] = True
                    for az in arange(360, az_static + self.width_step, self.width_step):
                        self.windows_is_gap[z][az-360] = True
            # debugging:
            print(z)
            for az in arange(0, 360, self.width_step):
                print("g " if self.windows_is_gap[z][az] else "1 " , end = '')
            print('')

base_x = 0
base_y = 0
base_z = 1
shift = 0.05

input_file_name = "forest_0.txt"

def calc_dist(x, y, z):
    return math.sqrt((x-base_x)*(x-base_x) + (y-base_y)*(y-base_y)
        + (z-base_z)*(z-base_z))

def calc_zenith(x, y, z):
    z = math.acos((z-base_z)/calc_dist(x, y, z))
    return math.degrees(z) + shift

def calc_azimuth(x, y, z):
    az = math.atan((y-base_y)/(x-base_x))
    if y - base_y < 0:
        az += math.pi
    elif x - base_x < 0:
        az += 2 * math.pi
    return math.degrees(az) + shift

def gap_equation(x, *arguments):
    gap, G, histogram = arguments
    total = 0
    n_bins = len(histogram)
    for i in arange(n_bins):
        total += math.exp(-G * x * (i+0.5)/n_bins) * histogram[i]
    return total - gap

def calculate_lai(x, cos_zenith, histogram):
    total = 0
    n_bins = len(histogram)
    for i in arange(n_bins):
        total += x * cos_zenith * i/n_bins * histogram[i]
    return total


verify_resolution = False
def filter(zenith, azimuth):
    if verify_resolution:
        if int(zenith * 10300501) % 10000 == 0:
            print(zenith)
        if zenith < 100.05 or zenith > 100.15:
            return False
    elif zenith < 0 or zenith > 90:
        return False
    return True
    

slice_height = 10
slice_step = 10
G = 0.5

def main():
    #f = open("PointsTREEElli12_01_1.txt", "r")
    f = open(input_file_name, "r")
    lines = f.readlines()
    f.close()

    windows = Windows()

    print('number of lines: {}\n'.format(len(lines)))
    azimuths = []
    for line in lines:
        numbers = line.split()
        x = float(numbers[0])
        y = float(numbers[2])
        z = float(numbers[1])

        zenith = calc_zenith(x, y, z)
        azimuth = calc_azimuth(x, y, z)
        if not filter(zenith, azimuth):
            continue
        windows.add_to_window(zenith, azimuth)
        azimuths.append(azimuth)
    if verify_resolution:
        azimuths = sorted(azimuths)
        for i in arange(len(azimuths)):
            print(azimuths[i])
    
    windows.tag_gap()
    total_lai = 0
    total_sin_z = 0
    for z in arange(0, 90, slice_step):
        start_z, end_z = windows.get_slice_boundary(z, slice_height)
        gap, histogram, big_gap_p = windows.get_gap_histogram(start_z, end_z)
        if gap == 1:
            continue
        print("zenith: {}, gap: {}".format(z, gap))
        x_guess = 1
        x = fsolve(gap_equation, x_guess, args=(gap, G, histogram))
        # print("rou * l_max: {}".format(x))
        cos_zenith = math.cos(math.radians(z + slice_step/2))
        sin_zenith = math.sin(math.radians(z + slice_step/2))
        lai = calculate_lai(x, cos_zenith, histogram)
        print("lai: {}, big_gap_p: {}".format(lai, big_gap_p))
        total_lai += lai * sin_zenith * big_gap_p
        total_sin_z += sin_zenith
    lai = total_lai / total_sin_z
    print("lai: {}".format(lai))

main()