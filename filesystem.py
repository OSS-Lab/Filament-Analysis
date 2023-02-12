import os
from pathlib import Path
import re

class DataFolders:
    def __init__(self):
        """ `DataFolders` intended use is as a convenient way to navigate data
        from the TIRF Experiments. Folder structure needs to be:
        (data_root)/Date/Sample_Name/Slides/Slide_Name/Files/
        Specify `data_root` in `working_dir/config.txt`
        """
        
        #: dict: {'Sample_1':  {'Slide_1': [['filepath_1', file_id_1], ...]}, ...}
        self.tree = {}
        #: str: We only want .tif files 
        self.default_extension = '.tif'
        #: path: config.txt --> Path(data_root). 
        self.base_dir = self._init_base_dir()
        #: list: List of Path with all dirs in data_root, where name ~ dd-mm-yyyy
        self.sub_folders = self._list_sub_folders()
        # str: dd-mm-yyyy string for current working dir
        self.date = self.sub_folders[0].name
        
        
    def _init_base_dir(self):
        """ Looks for config.txt and returns a Path object if found, none otherwise.  
        """
        current = Path(os.getcwd())
        config = current.joinpath('config.txt')
        
        if config.is_file():
            f = open(config, 'r')
            base_path = Path(f.readlines()[0])
            f.close()
        return base_path
        
        print('Initialisation Failed. Config.txt not found.')
        return None
    
    def _list_sub_folders(self):
        """ All folders with `dd-mm-yyyy` name type, within `base_dir`.
        """
        _, folders = self.qdir(self.base_dir)
        for item in folders:
            if re.search('\d\d-\d\d-\d{4}', str(item)) == None:
                folders.remove(item)
        return folders           
    
    @property
    def date(self):
        """ Get current date in dd-mm-yyyy format
        """
        return self._date
    
    @date.setter
    def date(self, date_string):
        """ Set current date from input string,in dd-mm-yyyy format
        """
        date_folder = self.base_dir.joinpath(date_string)
        if date_folder in self.sub_folders:
            self._date = date_folder.name
            self.tree = self.query_tree()
        else:
            print('Selected date is not a directory within', f'{self.base_dir}')

    @property
    def date_dir(self):
        """ Returns a pathlib.Path obj representing the current working date directory
        """
        return self.base_dir.joinpath(self.date)
       
    def qdir(self, dir_path):
        """Query a directoty. Needs a pathlib.Path as input
        
        Parameters
        ----------
        dir_path: path
            Current folder Path
        
        Returns
        -------
        files: list
            list of file Paths 
        folders: list
            list of sub-folder Paths
        """
        extension = self.default_extension
        files = []
        folders = []
        
        for item in dir_path.iterdir():
            if item.is_dir():
                folders.append(item)
            elif item.is_file():
                files.append(item)
        
        if extension!= '':
            files = [item for item in files if item.suffix == extension]
        
        return files, folders
    

    def query_tree(self):
        """Query the directory tree starting with the current `date`
        Returns
        -------
        dict:
            {'Sample_1':  {'Slide_1': [['filepath_1', file_id_1], ...]}, ...}
        """
        samples_dict = {}
        n_0 = 0
        
        _, sample_dirs = self.qdir(self.date_dir)
        
        for sample in sample_dirs:
            _, slide_dirs = self.qdir(sample.joinpath('Slides'))
            for slide in slide_dirs:
                files, _ = self.qdir(slide)
                files = [[file, n_0 + k] for k, file in enumerate(files)]
                n_0 += len(files)
                
                # ----- #
                if sample.name in samples_dict:
                    samples_dict[sample.name][slide.name] = files
                else:
                    samples_dict[sample.name] = {slide.name : files}
        return samples_dict

         
    def print_tree(self):
        """ Print the contents of `self.tree` in a more readable form.
        """
        if len(self.tree) == 0:
            print('Directory tree is empty, set a `date` folder first')
            return
        
        for sample in self.tree:
            print(sample)
            for slide in self.tree[sample]:
                print('\t', slide)
                for file in self.tree[sample][slide]:
                    print('\t','\t', file[0].name,f'  file_id ={file[1]}')
                    
    
    def path_from_id(self, id_num):
        """Returns the abs. path of a .tif stack, for a given id_num (int)
        """
        if len(self.tree) == 0:
            print('Directory tree is empty, set a `date` folder first')
            return
        # Naive search through for id, fix this another time
        for sample in self.tree:
            for slide in self.tree[sample]:
                files = self.tree[sample][slide]
                for file in files:
                    if file[1] == id_num:
                        return file[0]
        
        print('File ID Not Found')
        return None
                

    def axis_and_tracks(self, id_num):
        """Returns the abs. path of the filament axis and tracking results, for a given id_num (int)
        """
        datapath = str(self.path_from_id(id_num))
        datapath = re.sub(self.default_extension, '', datapath)

        axespath = re.sub('Slides', 'FilamentAxes', datapath)
        trackingpath = re.sub('Slides', 'Tracking', datapath)
        
        return Path(axespath), Path(trackingpath)
    
    def html_tree(self):
        """ Print the contents of `self.tree` with html tags.
        """
        if len(self.tree) == 0:
            print('Directory tree is empty, set a `date` folder first')
            return
        
        print('<ul>')
        
        for sample in self.tree:
            print('\t<li>')
            print('\t\t<details>')
            print(f'\t\t<summary> {sample} </summary>')
            
            print('\t\t<ul>')
            for slide in self.tree[sample]:
                print('\t\t\t<li>')
                print('\t\t\t<details>')
                # Directory name
                print('\t\t\t<summary>', slide, '</summary>')
                # Single dir file list
                print('\t\t\t<ul>')
                for file in self.tree[sample][slide]:
                    print('\t\t\t\t<li>', file[0].name, '</li>')
                print('\t\t\t</ul>')
                print('\t\t\t</details>')
                print('\t\t\t</li>')
            print('\t\t</ul>')
            print('\t\t</details></li>')
            
        print('</ul>')
    