import numpy as np
from scipy.spatial.distance import euclidean
from scipy.stats import linregress
from . import various_funcs as var


class SplineLine:
    
    def __init__(self, spline):
        """Base class for spline-like objects.

        Parameters
        ----------
        spline: numpy.ndarray 
        
        """
        self.spline = spline 
        self.hash_table = {}  # see compute_hash_table()
        self.search_fraction = 0.1 # see point_projection()
           
    def point_projection(self, point):
        """Finds the index of the spline point that best approximates an
        orthogonal projection of `point`.   
        
        Parameters
        ----------
        point : numpy.ndarray
            Point co-ordinates given as a vector of shape (1,2)
        
        Returns
        ------- 
        int
            Integer, index of the point on `self.spline` that best represents the projection of `point`
        """
        # Start at nearest index with regards to x-only
        nearest_ind = np.argmin(np.abs(self.spline[:,0] - point[0]))
        # Widen the search in the neighborhood. 
        search_radius = int(self.search_fraction*len(self.spline))
        # Account for proximity to edges
        lower = max(nearest_ind - search_radius, 0)
        upper = min(nearest_ind + search_radius, len(self.spline) - 1)
        
        distances = [euclidean(q,point) for q in self.spline[lower : upper,:]]
        
        return lower + np.argmin(distances)
    
    def xy_on_spline(self, point):
        """
        Project `point` onto `self.spline`. Calculates: 
        x := length of path (dist) from `self.spline`[0] 
        y := orthogonal distance from `point` to `self.spline`
        with +/- sign if `point` is above or below the line.
        
        Parameters
        ----------
        point: numpy.ndarray
            Point co-ordinates given as a vector of shape (1,2)
        
        Returns
        -------
        tuple
            Tuple of (array([x,y]) : numpy.ndarray, projection index : int)
        
        """
        # Index of best match spline projection
        proj_ind = self.point_projection(point)
        # Figure out if point above/below spline and adjust 'y'
        if point[1] > self.spline[proj_ind][1]:
            y = euclidean(point , self.spline[proj_ind])
        else:
            y = -euclidean(point , self.spline[proj_ind])
        # 'x' is just the distance from spline[0]
        if len(self.hash_table.keys()) !=0:
            x = self.hash_table[str(proj_ind)]
        else:
            raise ValueError("Call instance.compute_hash_table() first")
        
        return np.array([x,y]), proj_ind
    
    def compute_hash_table(self):
        """Compute a hash table (dict) to look up the distance (length of path)
        from spline[0] to any other point on the spline.
        """
        length = 0
        lengths_dict = {'0': 0}
        
        for k in range(1, len(self.spline)-1):
            length = length + euclidean(self.spline[k-1],self.spline[k])
            lengths_dict[str(k)] = length
        
        self.hash_table = lengths_dict
        
    
    def linear_ext(self, n, step, mode, to_hash = True):
        """
        Linear extrapolation of a spline end. Updates `self.spline`,
        and calls `self.compute_hash_table()` if to_hash = True.
        
        Parameters
        ----------
        n : int
            Number of pixels in the x-direction for which the spline will be extended.
        step : float
            dx for the points of the extension
        mode : str
            Use 'first' or 'last' to extend the leading or trailing end of spline
        to_hash : bool, optional
            Default True, pass false to supress hash, eg. for two or more extensions before
            you intend to calculate any projections or distances along the spline.
        """
        # Select appropriate pre/post x-range given `mode` and fit
        if mode == 'first':
            x = np.arange(self.spline[:,0][0] - n, self.spline[:,0][0],  step)
            a, b = linregress(self.spline[0:10,0], self.spline[0:10,1])[0:2]
        
        elif mode == 'last':
            x = np.arange(self.spline[:,0][-1], self.spline[:,0][-1] + n, step)
            a, b = linregress(self.spline[-10:,0], self.spline[-10:,1])[0:2]

        y = a*x + b
        ext = np.array([x,y]).T
        
        # Update self.spline and re-do the distance hash table
        if mode == 'first':
            self.spline = np.concatenate((ext, self.spline))
        else:
            self.spline = np.concatenate((self.spline, ext))
        
        if to_hash == True:
            self.compute_hash_table()



class FilamentAxis(SplineLine):
    
    def __init__(self, worm_axis, image_shape):
        """ `FilamentAxis` inehrits from `SplineLine`. 
        Made for working with axes of cyanobacterial filaments from TIRF.
        
        Parameters
        ----------
        image_shape : tuple
            Shape of original image
        """
        super().__init__(worm_axis)
        self.shape = image_shape
    
    def if_in_bounds(self):
        """Corrects spline if any computation puts it outside image bounds.
        """
        spl_x = self.spline[:,0]
        x_in_bounds = np.logical_and(spl_x >= 0, spl_x <= self.shape[1])
        if np.all(x_in_bounds) == True:
            return
        else:
            self.spline = self.spline[x_in_bounds,:]
        
        
    def extend(self):
        """Linear extension of both ends of filament axis, 
        until the boundaries of the image.
        """
        # dx for extrapolated ends
        step = np.abs(np.round(self.spline[1][0] - self.spline[0][0],4))
        
        self.if_in_bounds()
        
        # number of pixels in 'x' direction 
        n_pre = int(self.spline[0][0]) 
        n_post = int(self.shape[1] - self.spline[-1][0])
        # Extend both ends, without computing hash table in between
        self.linear_ext(n_pre, step, \
                        mode = 'first', to_hash = False)
        
        self.linear_ext(n_post, step, \
                        mode = 'last', to_hash = False)
        # Compute hash table
        self.compute_hash_table()
    
    def shift(self, distance):
        """Perpendicularly shift every point of the spline by `distance`"""
        all_shifted = []
        
        for k in range(0,len(self.spline)):
            if k !=0 and k!=len(self.spline)-1:
                tangent = self.spline[k+1] - self.spline[k-1]
            elif k==0:
                tangent = self.spline[1] - self.spline[0]
            else:
                tangent = self.spline[-1]-self.spline[-2]
                
            point_shifted = self.spline[k] + var.perpendicular(tangent, np.abs(distance), np.sign(distance))
            all_shifted.append(point_shifted)
        
        self.spline = np.array(all_shifted)
        self.if_in_bounds()
        
    @classmethod
    def from_image(cls, input_image, shift_val = 0, prune_n = 7, spline_k = 4):
        """Alternative constructor that works with a binary input
        image rather than a spline array.

        Parameters
        ----------
        input_image : numpy.ndarray 
        shift_val : int
            Default 0, use if the spline needs shifting.
        prune_n : int
            How many pixels should the branches of the image skeleton be pruned
            for. Default 7.
        spline_k : int
            Degree of the smoothing spline. Default 4.
        """
        spline = var.filament_spline(input_image, prune_n, spline_k)
        spline_obj = cls(spline, input_image.shape)
        if shift_val != 0:
            spline_obj.shift(shift_val)
        spline_obj.extend()
        return spline_obj
