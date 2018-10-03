typeset -i k=0
INICIO=81
FINAL=170
for ((k=$INICIO;k<=$FINAL;k++))
#for k in {1..3}
do
mv -f /Disco2/Pheno/Backgrounds/W+jets/muestra81-160/W+jets_$k/ /Disco2/Pheno/Backgrounds/W+jets/  #Change folder name.
#rm -r /Disco2/Pheno/Backgrounds/W+jets/W+jets_$k/ 
echo File W+jets_$k moved.       #Change folder name.
done

