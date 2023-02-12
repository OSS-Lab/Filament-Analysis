import numpy as np
from Image_Analysis.various_funcs import overlap_len

class PandasTracks:
    
    def __init__(self, dataframe):
        """Interfaces with a `pandas.DataFrame` to conviniently extract
        data and perform operations on it.
        
        Parameters
        ----------
        dataframe : pandas.DataFrame
        
        """
        self.data = dataframe
        self.track_ids = self.all_track_ids()

    def all_track_ids(self):
        """ Fetch all unique track ID's from the data frame.
        Stores them as `self.track_ids`
        """
        track_ids = np.unique(self.data['TRACK_ID'].to_numpy(dtype=np.int16))   
        track_beginnings = np.array([])
        
        for tr_id in track_ids:
            str_id = str(tr_id)
            frames = self.get_frames(str_id)
            track_beginnings = np.append(track_beginnings,frames[0])
        
        sorted_args = np.argsort(track_beginnings)
        track_ids = track_ids[sorted_args]
        
        return [str(x) for x in track_ids]
    
    
    def frame_spots(self, frame_id):
        """ Fetch all spots present within a given frame
        
        Parameters
        ----------
        frame_id : int
        
        Returns
        -------
        pandas.DataFrame
            (ID Column) | POSITION_X | POSITION_Y | TRACK_ID 
          
        """
        
        
        all_spots_in_frame = self.data.loc[self.data['FRAME'] == str(frame_id),\
                                           ['POSITION_X', 'POSITION_Y', 'TRACK_ID']]
        
        return all_spots_in_frame
    
    def get_tracks(self, tr_id_list):
        """ For a list of track ID's, get all the points corresponding to these tracks
        
        Parameters
        ----------
        tr_id_list : list of str
        
        Returns
        -------
        dict
            Dictionary {str: numpy.ndarray} where keys are track_id's and the array 
            columns are X, Y, Frame
        """
        output = {}
        for tr_id in tr_id_list:
            x_y_frame = self.data.loc[self.data['TRACK_ID'] == tr_id, \
                                      ['POSITION_X', 'POSITION_Y', 'FRAME']]
            
            x_y_frame = x_y_frame.to_numpy(dtype = np.float64)
            sorted_time_args = np.argsort(x_y_frame[:,-1])
            x_y_frame = x_y_frame[sorted_time_args]
            output[tr_id] = x_y_frame
        return output
    
    def get_frames(self, tr_id):
        """ Get all frames in sorted order for a given track ID 
        Parameters
        ----------
        tr_id : str
        
        Returns
        -------
        numpy.ndarray
        """
        frames = self.data.loc[self.data['TRACK_ID'] == tr_id, \
                               ['FRAME']].to_numpy(dtype=np.int16)
        frames = frames.reshape(len(frames),)
        
        return np.sort(frames)
        
    def tracks_endpoints(self):
        """ Get Min Frame and Max Frame for every trajectory ID.
        Returns
        -------
        dict
            {tr_id : [min, max]}
        """
        min_max = {}
        for tr_id in self.track_ids:
            frames = self.get_frames(tr_id)
            min_max[tr_id] = [frames[0], frames[-1]]
        return min_max

    
    def overlap_matrix(self, overlap_threshold):
        """Compute overlap matrix """
        min_max = {}
        for tr_id in self.track_ids:
            frames = self.get_frames(tr_id)
            min_max[tr_id] = [frames[0], frames[-1]]
        
        size = len(self.track_ids)
        ovp_matrix = np.zeros((size,size)).astype(float)
        
        for n,p in enumerate(self.track_ids):
            for m,q in enumerate(self.track_ids):
                overlap_value = overlap_len(min_max[p],min_max[q])
                if overlap_value != None:
                    if overlap_value > overlap_threshold:
                        ovp_matrix[n , m] = overlap_value

        return (ovp_matrix > overlap_threshold).astype(int)