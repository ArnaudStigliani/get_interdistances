import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import pprint
from BCBio.GFF import GFFExaminer

in_file = "../../extraPrograms_or_big_files/TAIR10_GFF3_genes.gff"
examiner = GFFExaminer()
in_handle = open(in_file)
pprint.pprint(examiner.parent_child_map(in_handle))
in_handle.close()