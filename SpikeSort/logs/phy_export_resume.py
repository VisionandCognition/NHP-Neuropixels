import shutil, spikeinterface.full as si
from spikeinterface.exporters import export_to_phy
S="/media/DATA1/spikesort_scratch/m032_20251010_exp2/"
O="/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_analysis/SpikeSorting/m032/20251010/exp2/phy"
an=si.load_sorting_analyzer(S+"analyzer")
export_to_phy(an,S+"phy",compute_pc_features=False,copy_binary=True,remove_if_exists=True,verbose=False)
print("exported locally, copying")
shutil.copytree(S+"phy",O,dirs_exist_ok=True)
print("copied")
