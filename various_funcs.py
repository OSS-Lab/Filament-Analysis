import sys
import numpy as np
sys.path.insert(0, 'C:\\Warwick\\Code\\')
import diocane_1 as dio
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize
from scipy import interpolate
from skimage.measure import label, regionprops
from scipy.signal import convolve2d



def filament_spline(img, prune_n = 7, spline_k = 4):

    # Construct a skeleton
    kernels = dio.load_kernels()
    skl = skeletonize(img.astype(np.uint8), method = 'lee')
    pruned = dio.prune_skeleton(skl, kernels, prune_n)

    # Get the (X,Y) points of the pruned skeleton
    data = np.argwhere(pruned[0]>0)

    # Split X and Y
    spl_x = data[:,1]
    spl_y = data[:,0]

    # Sort indices by X value
    sorted_ind = np.argsort(spl_x)

    # Apply the sorted index
    spl_x = spl_x[sorted_ind]
    spl_y = spl_y[sorted_ind]

    # Fit spline to this skeleton
    spline = interpolate.UnivariateSpline(spl_x,spl_y, k = spline_k)
    x_interp = np.arange(0, np.max(spl_x), 0.1)
    y_interp = spline(x_interp)
    
    spl_twocol = np.array([x_interp, y_interp]).T
    
    return spl_twocol

def simple_threshold(image, thr_val):
    binary = image > thr_val
    return binary.astype(np.uint8)


    
def filter_threshold(image, thr_val, min_blob_area):
    binary = image > thr_val
    binary =  binary.astype(np.uint8)
    labelled_binary = label(binary)
    props = regionprops(labelled_binary)
    to_remove =[p.label for p in props if p.area<min_blob_area]
    for item in to_remove:
        binary[labelled_binary == item] = 0
    return binary
    


def two_thresholds(image, thr_bottom, thr_top):
    binary_1 = image > thr_bottom 
    binary_2 = image < thr_top
    both = np.logical_and(binary_1, binary_2)
    return both.astype(np.uint8)

def image_histogram(image, bins=64, xlims=[], ylims=[]):
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    ax.hist(image.ravel(), bins = bins)
    if len(xlims)!=0:
        ax.set_xlim(xlims[0], xlims[1])
    if len(ylims)!=0:
        ax.set_ylim(ylims[0], ylims[1])
    plt.show()

def plot_two(image_1, image_2, spline_1 = [], spline_2 = []):
    fig, ax = plt.subplots(1,2, figsize=(12,8))
    ax[0].imshow(image_1, cmap='gray')
    ax[1].imshow(image_2, cmap='gray')
    if len(spline_1) != 0:
        ax[0].plot(spline_1[:,0], spline_1[:,1], '-', color='red', linewidth=0.4)
    if len(spline_2) != 0:
        ax[1].plot(spline_2[:,0], spline_2[:,1], '-', color='yellow', linewidth=0.4)
    plt.show()
    
def exclude(prev_filament, image):
    to_exclude = np.argwhere(prev_filament > 0)
    excluded = np.copy(image)
    for point in to_exclude:
        excluded[point[0], point[1]] = 0
    return excluded

def spline_erase(spline, image, mode):
    image_copy = np.copy(image)
    spline_copy = np.copy(spline).astype(int)
    uniques = []
    # There will be duplicates in int(spline)
    for k in range(0, len(spline_copy)):
        point  =  spline_copy[k].tolist()
        if point not in uniques:
            uniques.append(point)
    # some splines may run past image borders, i think i fixed this in filaxis cls
    # edit later
    #for k in range(0, len(uniques)):
        #if uniques[k][0] < image_copy.shape[1]:
            #if uniques[k][1] < image_copy.shape[0]:
            
    
    for point in uniques:
        if mode == 'over':  
            to_del = range(point[1],image_copy.shape[0])                       
        elif mode == 'under': 
            to_del = range(point[1], -1, -1)
        for n in to_del:
            image_copy[n, point[0]] = 0
    
    return image_copy

def perpendicular(vector, distance, orientation):
    rot = orientation*np.array([[0,-1],[1,0]]) # +90 rot
    
    normal = np.dot(rot, vector.reshape(2,1))
    normal = normal/np.linalg.norm(normal)
    
    return (distance*normal).reshape(2,)

def show_image(image, splines = []):
    fig, ax = plt.subplots(1,1, figsize = (10,8))
    colors = ['red', 'cyan', 'orange', 'pink', 'purple']
    ax.imshow(image, cmap='gray')
    if len(splines) > 0:
        for k, spline in enumerate(splines):
            ax.plot(spline.spline[:,0], spline.spline[:,1], '-', color=colors[k], linewidth=1)
    plt.show()
        
# this below needs cleaning




def get_endpoints(differences_dict):
    traj_endpoints = {}
    for track_id in differences_dict.keys():
        mid_points = []
        for item in differences_dict[track_id]:
            mid_points.append(item[-2] + (item[-1] - item[-2])/2)
        traj_endpoints[track_id] = [min(mid_points), max(mid_points)]
    return traj_endpoints
    

def overlap_matrix(traj_endpoints,overlap_threshold):
    traj_ids = traj_endpoints.keys()
    size = len(traj_ids)
    ovp_matrix = np.zeros((size,size)).astype(float)
    
    for n,p in enumerate(traj_ids):
        for m,q in enumerate(traj_ids):
            overlap_value = overlap_len(traj_endpoints[p],traj_endpoints[q])
            if overlap_value != None:
                if overlap_value > overlap_threshold:
                    ovp_matrix[n , m] = overlap_value
    
    return (ovp_matrix > overlap_threshold).astype(int)

def continous_counts(some_array):
    lengths = []
    counter = 0
    for value in some_array:
        if value == 1:
            counter = counter + 1
        else:
            if counter!=0:
                lengths.append(counter)
            counter = 0
    if counter !=0:
        lengths.append(counter)
    return lengths
        
        
    
def largest_continuum(ovp_matrix):
    size = len(ovp_matrix)
    row_lengths = []
    col_lengths = []
    for n in range(0, size):
        row_counts = continous_counts(ovp_matrix[n, :])
        col_counts = continous_counts(ovp_matrix[:, n])
        row_lengths += row_counts
        col_lengths += col_counts
    return np.min([np.max(row_lengths), np.max(col_lengths)])
        
    
def sub_matrices(ovp_matrix):
    """needs bool """
    size_0 = largest_continuum(ovp_matrix)
    sizes = [n for n in range(1, size_0 + 1)]
    candidates = {}
    
    for size in sizes:
        kernel = np.ones((size, size)).astype(int)
        convolved = (1/size**2)*convolve2d(ovp_matrix, kernel, mode = 'valid')
        convolved_bool = convolved == 1
        if np.any(convolved_bool):
            max_args = np.argwhere(convolved_bool == True)
            for n in range(0, len(max_args)):
                row, col = max_args[n,:]
                key_size = str(size)
                if key_size in candidates.keys():
                    candidates[key_size].append([row, col])
                else:
                    candidates[key_size] = [[row, col]]
        
    sizes = sorted([int(n) for n in candidates.keys()])
    
    if len(sizes) > 3:
        sizes = sizes[-3:]
    else:
        sizes = sizes[-1:]
    
    sizes.reverse()
    
    offs_all = []
    
    for size in sizes:
        rows_cols = candidates[str(size)]
        offsets = [item[0] for item in rows_cols if item[1]==item[0]]
        good_offsets = [offsets[0]]
        n = 0
        
        for offset in offsets:
            if offset - good_offsets[n] > size:
                good_offsets.append(offset)
                n = n + 1
        
        offs_all.append(good_offsets)
    
    if len(sizes) > 1:
        for n in range(1, len(sizes)):
            all_prev_offs = [q for p in offs_all[:n] for q in p]
            offs_all[n] = [item for item in offs_all[n] if item not in all_prev_offs]
    
    fully_processed = {}
    
    for k, size in enumerate(sizes):
        if len(offs_all[k]) > 0:
            key = str(size)
            values = offs_all[k]
            fully_processed[key] = values
    
    return fully_processed


def overlapping_tracks(track_ids, ovp_matrix):
    sections = sub_matrices(ovp_matrix)
    tracks = list(track_ids)
    groups = []
    for size in sections.keys():
        for offset in sections[size]:
            groups.append(tracks[offset : offset + int(size)])
    return groups



def check_num(x):
    if type(x[0]) == type('abc'):
        return [int(x[0]), int(x[1])]
    else:
        return x
#
def overlap_len(a,b):
    a, b = check_num(a), check_num(b)
    d_left, d_right = a[0] - b[1], a[1] - b[0]
    a_len, b_len = a[1] - a[0], b[1] - b[0]
    
    if d_left >= 0 or d_right<=0:
        return None
    
    if d_left < 0 and -d_left <= a_len:
        if a[0] >= b[0]:
            return round((b[1] - a[0])/a_len,2)
        else:
            return round(b_len / a_len,2)
    
    if d_right > 0 and d_right <= a_len:
        if b[1] > a[1]:
            # match with right overhang
            return round((a[1] - b[0])/a_len,2)

    elif d_right > 0 and d_right > a_len:
        #b encompasses a
        if -d_left > a_len:
            return 1.0
        
    # If no cases were found something went off
    print('Something went wrong')
    return None

def perpendicular(vector, distance, orientation):
    rot = orientation*np.array([[0,-1],[1,0]]) # +90 rot
    
    normal = np.dot(rot, vector.reshape(2,1))
    normal = normal/np.linalg.norm(normal)
    
    return (distance*normal).reshape(2,)
    
def match_x(x_0, x_array):
    x_diff = np.abs(x_array - x_0)
    arg = np.argmin(x_diff)
    return int(arg)
    
    


def overlap_matrix_2(frames_list, overlap_threshold):
    """ Simpler overlap matrix for final processing"""
    min_max = []
        
    for frames in frames_list:
        min_max.append([min(frames), max(frames)])
 
    size = len(frames_list)
    ovp_matrix = np.zeros((size,size)).astype(float)
        
    for i in range(len(min_max)):
        for j in range(len(min_max)):
    
            overlap_value = overlap_len(min_max[i],min_max[j])
            if overlap_value != None:
                if overlap_value > overlap_threshold:
                    ovp_matrix[i,j] = overlap_value

    return (ovp_matrix > overlap_threshold).astype(int)
  
    