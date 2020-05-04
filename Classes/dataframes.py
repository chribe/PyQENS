from Measurement import Data
import mantid.simpleapi as MTD
import numpy as np
import pickle

filepath_raw='\\\\moneta.uni-tuebingen.de\\groupshare\\Biophysics\\group_share\\Beamtime\\ILL-IN16B-9-13-829-22022020\\exp_9-13-829\\rawdata\\'
filepath_analysed='\\\\moneta.uni-tuebingen.de\\groupshare\\Biophysics\\group_share\\Beamtime\\ILL-IN16B-9-13-829-22022020\\exp_9-13-829\\processed\\qens_analysed\\'
filepath_van='\\\\moneta.uni-tuebingen.de\\groupshare\\Biophysics\\group_share\\Beamtime\\ILL-IN16B-9-13-829-22022020\\exp_9-13-829\\rawdata_829-879\\'

MTD.config.appendDataSearchDir( filepath_raw )
MTD.config.appendDataSearchDir( filepath_analysed )
MTD.config.appendDataSearchDir( '\\\\moneta.uni-tuebingen.de\\groupshare\\Biophysics\\group_share\\Beamtime\\ILL-IN16B-9-13-829-22022020\\exp_9-13-829\\rawdata\\rawdata\\')

vanadium=''

for r in range(268794,268799):
    vanadium = vanadium+filepath_van+str(r)+'.nxs,'

vanadium=vanadium[:-1]


EC=''

for r in range(265619,265627):
    EC = EC+filepath_van+str(r)+'.nxs,'

EC=EC[:-1]

#to load data use: Data.load_data(pathraw, filepath_analysed, measurement, Vanadium, EC)
Beamtime=[] # create an open list to put in the files

Beamtime += [Data(namesample='PEG2.5%', run_num=list(range(288872, 288882)), type=['qens'])]#291
Beamtime += [Data(namesample='PEG2.5%', run_num=list(range(288834, 288842)), type=['qens'])]#280
Beamtime += [Data(namesample='PEG5%',   run_num=list(range(288947, 288957)), type=['qens'])]#291

Beamtime += [Data(namesample='IgPEG6%21Cdense',   run_num=list(range(290364, 290372)), type=['qens'])]#277
Beamtime += [Data(namesample='IgPEG6%21Cdense',   run_num=list(range(291897, 291905)), type=['qens'])]#273
Beamtime += [Data(namesample='IgPEG6%21Cdense',   run_num=list(range(285895, 285919)), type=['qens'])]#291


Beamtime += [Data(namesample='IgPEG8%21Cdense',   run_num=list(range(292022, 292030)), type=['qens'])]#310
Beamtime += [Data(namesample='IgPEG8%21Cdense',   run_num=list(range(292635, 292643)), type=['qens'])]#280


Beamtime +=[Data(namesample='D2O',   run_num=list(range(268799, 268809)), type=['qens'])]#310

for meas in Beamtime:
    meas.load_data(filepath_van, filepath_analysed, Vanadium=vanadium, EC=EC)


