import math

base_x = 9
base_y = 15
base_z = 1

interval = 0.04

def calc_dist(x, y, z):
    return math.sqrt((x-base_x)*(x-base_x) + (y-base_y)*(y-base_y)
        + (z-base_z)*(z-base_z))

def main():    
    f = open("PointsTREEElli12_01_1.txt", "r")
    lines = f.readlines()
    nearest = 0
    max_dist = 0.0
    distances = []
    pathlen_bins = []

    for i in range(int(1 / interval) + 1):
        pathlen_bins.append(0)

    print('number of lines: {}\n'.format(len(lines)))
    for line in lines:
        numbers = line.split()
        if (numbers[-1] == '1'):
            continue
        if (numbers[-2] == '1'):
            nearest = calc_dist(float(numbers[0]), float(numbers[1]), float(numbers[2]))
        elif (numbers[-2] == numbers[-1]):
            if nearest == 0:
                # print(numbers[-1])
                continue
            farthest = calc_dist(float(numbers[0]), float(numbers[1]), float(numbers[2]))
            distances.append(farthest - nearest)
            if farthest - nearest > max_dist:
                max_dist = farthest - nearest
            nearest = 0
    for i, d in enumerate(distances):
        relative_d = d / max_dist
        index = int(relative_d / interval)
        pathlen_bins[index] += 1

    for d in pathlen_bins:
        print(d)


def calc_zenith(x, y, z):
    return math.acos((z-base_z)/calc_dist(x, y, z))

def calc_azimuth(x, y, z):
    return math.atan((y-base_y)/(x-base_x))

window_width = 0.04
window_height = 0.04
full_count = 4 * window_width * window_height / (0.0016 * 0.0016)
step = 0.005
count_threshold = 2

window_width_p = 0.02
window_height_p = 0.02
step_p = 0.005
ignore_pathlen_above = 25
weight_of_adjustment = 0.5
edge_ratio = 0.1


def with_all_rays():
    f = open("PointsTREEElli12_01_1.txt", "r")
    lines = f.readlines()
    f.close()

    # add multi echo rays
    f = open("PointsTREEElli12_01_m.txt", "r")
    other_lines = f.readlines()
    for line in other_lines:
        numbers = line.split()
        if (numbers[-2] == numbers[-1]):
            lines.append(line)

    max_pathlen = 0.0
    max_count = 0
    zenith_azimuth = {}
    zenith_dist = {}
    zenith_pathlen = {}
    pathlens = []
    pathlen_bins = []
    cnt_bins = []
    neg_gap_p_bins = []
    final_bins = []
    zeniths = []
    azimuths = []
    zenith_min = 100
    azimuth_min = 100
    zenith_max = -100
    azimuth_max = -100

    for i in range(int(1 / interval) + 1):
        pathlen_bins.append(0)
    for i in range(int(1 / interval) + 1):
        neg_gap_p_bins.append(0)
    for i in range(int(1 / interval) + 1):
        cnt_bins.append(0)
    for i in range(int(1 / interval) + 1):
        final_bins.append(0)

    print('number of lines: {}\n'.format(len(lines)))
    for line in lines:
        numbers = line.split()
        x = float(numbers[0])
        y = float(numbers[1])
        z = float(numbers[2])
        zenith = calc_zenith(x, y, z)
        azimuth = calc_azimuth(x, y, z)
        zeniths.append(zenith)
        azimuths.append(azimuth)
        if zenith < zenith_min:
            zenith_min = zenith
        if zenith > zenith_max:
            zenith_max = zenith
        if azimuth < azimuth_min:
            azimuth_min = azimuth
        if azimuth > azimuth_max:
            azimuth_max = azimuth

        if zenith in zenith_azimuth:
            print("warning: replicated key")
        else:
            zenith_azimuth[zenith] = azimuth
        if zenith in zenith_dist:
            print("warning: replicated key")
        else:
            zenith_dist[zenith] = calc_dist(x, y, z)
    zeniths = sorted(zeniths)
    azimuths = sorted(azimuths)
    zenith_low = zeniths[int(len(zeniths)*edge_ratio)]
    zenith_high = zeniths[int(len(zeniths)*(1-edge_ratio))]
    azimuth_low = azimuths[int(len(azimuths)*edge_ratio)]
    azimuth_high = azimuths[int(len(azimuths)*(1-edge_ratio))]
    print(zenith_low)
    print(zenith_high)
    print(azimuth_low)
    print(azimuth_high)

    y = zenith_min
    y_index = 0
    while y < zenith_max:
        x = azimuth_min
        while x < azimuth_max:
            # start a window
            y0_index = y_index
            dist_min = -1
            dist_max = -1
            count = 0
            while y0_index < len(zeniths) and zeniths[y0_index] - y < window_height:
                x0 = zenith_azimuth[zeniths[y0_index]]
                if x0 > x and x0 - x < window_width:
                    dist = zenith_dist[zeniths[y0_index]]
                    if dist < dist_min or dist_min == -1:
                        dist_min = dist
                    if dist > dist_max:
                        dist_max = dist
                    count += 1
                y0_index += 1
            # finish a window
            if count >= count_threshold:
                pathlen = dist_max - dist_min
                if max_pathlen < pathlen:
                    max_pathlen = pathlen
                pathlens.append(pathlen)
                if count >= max_count:
                    max_count = count
                # rerun the window
                y0_index = y_index
                while y0_index < len(zeniths) and zeniths[y0_index] - y < window_height:
                    x0 = zenith_azimuth[zeniths[y0_index]]
                    if x0 > x and x0 - x < window_width:
                        zenith_pathlen[zeniths[y0_index]] = pathlen
                    y0_index += 1
            x += step
        y += step
        while y_index < len(zeniths) and zeniths[y_index] < y:
            y_index += 1

    print("max count: {}, full count: {}".format(max_count, full_count))
    for i, d in enumerate(pathlens):
        relative_d = d / max_pathlen
        index = int(relative_d / interval)
        pathlen_bins[index] += 1

    print("path lens:")
    for d in pathlen_bins:
        print(d)
    
    y = zenith_min
    y_index = 0
    while y < zenith_max:
        x = azimuth_min
        while x < azimuth_max:
            # start a window
            y0_index = y_index
            count = 0
            avg_pathlen = 0
            if x < azimuth_low or x + window_width_p > azimuth_high\
                    or y < zenith_low or y + window_height_p > zenith_high:
                x += step_p
                continue
            while y0_index < len(zeniths) and zeniths[y0_index] - y < window_height_p:
                x0 = zenith_azimuth[zeniths[y0_index]]
                if x0 > x and x0 - x < window_width_p:
                    if not zeniths[y0_index] in zenith_pathlen:
                        y0_index += 1
                        continue
                    count += 1
                    avg_pathlen += zenith_pathlen[zeniths[y0_index]]
                y0_index += 1
            # finish a window
            if count >= 1:
                neg_gap_p = -math.log(1 - count / full_count)
                avg_pathlen = avg_pathlen / max_pathlen / count / interval
                neg_gap_p_bins[int(avg_pathlen)] += neg_gap_p
                cnt_bins[int(avg_pathlen)] += 1
            x += step_p
        y += step_p
        while y_index < len(zeniths) and zeniths[y_index] < y:
            y_index += 1
    
    for i in range(len(neg_gap_p_bins)):
        if cnt_bins[i] != 0:
            neg_gap_p_bins[i] /= cnt_bins[i]

    
    y = zenith_min
    y_index = 0
    finals = []
    max_final = 0
    while y < zenith_max:
        x = azimuth_min
        while x < azimuth_max:
            # start a window
            y0_index = y_index
            count = 0
            avg_pathlen = 0
            is_edge = False
            if x < azimuth_low or x + window_width_p > azimuth_high\
                    or y < zenith_low or y + window_height_p > zenith_high:
                is_edge = True
            while y0_index < len(zeniths) and zeniths[y0_index] - y < window_height_p:
                x0 = zenith_azimuth[zeniths[y0_index]]
                if x0 > x and x0 - x < window_width_p:
                    if not zeniths[y0_index] in zenith_pathlen:
                        y0_index += 1
                        continue
                    count += 1
                    avg_pathlen += zenith_pathlen[zeniths[y0_index]]
                y0_index += 1
            # finish a window
            if count >= 1:
                neg_gap_p = -math.log(1 - count / full_count)
                avg_pathlen = avg_pathlen / max_pathlen / count / interval
                avg_neg_gap_p = neg_gap_p_bins[int(avg_pathlen)]
                if is_edge or int(avg_pathlen) > ignore_pathlen_above:
                    final = avg_pathlen
                else:
                    # print(neg_gap_p / avg_neg_gap_p)
                    final = avg_pathlen * neg_gap_p / avg_neg_gap_p * weight_of_adjustment\
                            + avg_pathlen * (1-weight_of_adjustment)
                if final > max_final:
                    max_final = final
                finals.append(final)
            x += step_p
        y += step_p
        while y_index < len(zeniths) and zeniths[y_index] < y:
            y_index += 1
    for i, d in enumerate(finals):
        relative_d = d / max_final
        index = int(relative_d / interval)
        final_bins[index] += 1

    print("final pathlen distribution:")
    for d in final_bins:
        print(d)

with_all_rays()         
