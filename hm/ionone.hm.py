
from modeller import *
from modeller.automodel import *

# log.verbose()
env = Environ()

class MyModel(DOPEHRLoopModel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
        rsr.add(secondary_structure.Alpha(self.residue_range('18:A', '49:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('57:A', '84:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('92:A', '127:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('137:A', '162:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('191:A', '227:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('230:A', '262:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('268:A', '294:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('296:A', '320:A')))
        rsr.add(secondary_structure.Alpha(self.residue_range('180:A', '186:A')))

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
env.io.hetatm = True

a = DOPEHRLoopModel(env,
              alnfile  = 'ionone.ali',
              knowns   = ('OR51E2.beta-ionone'),
              sequence = 'OR51E2')
a.starting_model = 0
a.ending_model   = 9
a.library_schedule = autosched.slow
a.max_var_iterations = 300
# a.md_level = refine.very_slow
a.make()
