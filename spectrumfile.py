import re
from pyteomics import mzml, mgf
from os.path import splitext
from collections import defaultdict

class SpectrumFile:
    def __init__(self, file_path):
        self.indices = defaultdict(dict)
        self.spectra_source = None
        self.file_format = None
        self._build_index = None
        self.get_by_id = None
        
        self._load(file_path)
        self._build_index()
        
    def _load(self, file_path):
        extension = splitext(file_path)[1]

        if extension.lower() == ".mzml":
            print(f"Inferred mzML format from {file_path}")
            self.spectra_source = mzml.MzML(file_path)
            self.file_format = 'mzml'
            self._build_index = self._index_MZML
            self.get_by_id = self._get_by_id_MZML
            
        elif extension.lower() == ".mgf":
            print(f"Inferred MGF format from {file_path}")
            self.spectra_source = mgf.IndexedMGF(file_path)
            self.file_format = 'mgf'
            self._build_index = self._index_MGF
            self.get_by_id = self._get_by_id_MGF
        
        else:
            print(
                f"Cannot infer format from {file_path}, only mzML and MGF formats are supported"
            )
            raise Exception("Unsupported spectra file format")
    
    def _get_by_id_MGF(self, id_string):
        if type(id_string) is int:
            return self.spectra_source.get_by_index(self.indices['scan'][id_string])
        elif type(id_string) is str:
            #in case of whitespaces
            id_string = id_string.strip()
            try:
                #the default way
                return self.spectra_source.get_by_id(id_string)
            except KeyError:
                #trying to recover
                match_index = re.match(r'index(?:\s+)?=(?:\s+)?(\d+)', id_string)

                if not match_index is None:
                    return self.spectra_source.get_by_index(int(match_index.group(1)))

                if id_string in self.indices['title'].keys():
                    return self.spectra_source.get_by_index(self.indices['title'][id_string])

                if type(id_string) is int or id_string.isdigit():
                    return self.spectra_source.get_by_index(self.indices['scan'][int(id_string)])

                raise KeyError(f'Cannot infer MGF scan from spectrum_id={id_string}')
                
        else:
            raise TypeError(f'Unsupported id type ({type(id_string)}): should be int or string')

    def _get_by_id_MZML(self, id_string):
        if type(id_string) is int:
            return self.spectra_source.get_by_index(self.indices['scan'][id_string])
        elif type(id_string) is str:
            #in case of whitespaces
            id_string = id_string.strip()
            try:
                #the default way
                return self.spectra_source.get_by_id(id_string)
            except KeyError:
                #trying to recover
                match_index = re.match(r'index(?:\s+)?=(?:\s+)?(\d+)', id_string)

                if not match_index is None:
                    return self.spectra_source.get_by_index(int(match_index.group(1)))

                #if id_string in self.indices['title'].keys():
                #    return self.spectra_source.get_by_index(self.indices['title'][id_string])

                if type(id_string) is int or id_string.isdigit():
                    return self.spectra_source.get_by_index(self.indices['scan'][int(id_string)])

                raise KeyError(f'Cannot infer MzML scan from spectrum_id={id_string}')
                
        else:
            raise TypeError(f'Unsupported id type ({type(id_string)}): should be int or string')
    
    def _index_MGF(self):
        for index in range(len(self.spectra_source)):
            params = self.spectra_source[index]['params']
            if 'title' in params.keys():
                self.indices['title'][params['title']] = index
                
            if 'scans' in params.keys():
                self.indices['scan'][int(params['scans'])] = index
    
    def _index_MZML(self):
        for spectrum in self.spectra_source:
            spectrumID = spectrum['id']
            index = spectrum['index']
            
            scan_match = re.match(r'.+scan(?:\s+)?=(?:\s+)?(\d+)', spectrumID)
            if not scan_match is None: 
                self.indices['scan'][int(scan_match.group(1))] = index
            
            self.indices['id'][spectrumID] = index
       