typeset -i k=0
INICIO=0   
FINAL=300
for ((k=$INICIO;k<=$FINAL;k++))
#for k in {1..3}
do
# cd /Disco2/Pheno/Backgrounds/W+jets-hadronic/W+jets-hadronic_$k/Events/run_01   #Change
#cd /Disco1/Pheno/BackgroudSamples/multijets/multijets_$k/Events/run_01
#cd /Disco2/Pheno/Backgrounds/DY+jets/DY+jets_$k/Events/run_01
cd /Disco1/Pheno/BackgroudSamples/singletop/singletop_$k/Events/run_01
#cd /Disco1/Pheno/monotop/monotop_full/monotop_full_$k/Events/run_01
#cd /Disco1/Pheno/monotop/monotop1jet_full/monotop1jet_full_$k/Events/run_01
#cd /Disco1/Pheno/monotop/monotop1jet_left/monotop1jet_left_$k/Events/run_01
#cd /Disco1/Pheno/monotop/monotop1jet_right/monotop1jet_right_$k/Events/run_01
rm -r output_pythia8.root
rm -r output_pythia6.root 
rm -r events.lhe.gz
rm -r run_01_tag_1_banner.txt
rm -r tag_1_pythia_beforeveto.tree.gz
rm -r tag_1_pythia_events.lhe.gz
rm -r tag_1_pythia_events.tree.gz
rm -r tag_1_pythia_lhe_events.root
rm -r tag_1_pythia.log
rm -r tag_1_pythia.root
rm -r tag_1_pythia_xsecs.tree
rm -r unweighted_events.root
rm -r unweighted_events.lhe
rm -r tag_1_pythia_events.hep
echo Delete files in singletop_$k.       #Change folder name.
done
