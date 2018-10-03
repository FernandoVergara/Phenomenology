typeset -i i=0
for x in `cat names.dat` ; do
mv $x SingleMuPt10_pythia8_cfi_py_GEN_SIM_DIGI_$i.root
echo $x
echo $i
i=i+1
done

