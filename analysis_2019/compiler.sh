typeset -i i=0
typeset -i j=0
typeset -i k=0
typeset -i start=1
typeset -i end=1
typeset -i delta=10

	make compile_ROOT_Delphes
####################################################################################
	typeset -i startsignal1=1
        typeset -i endsignal1=1
        for ((j=$startsignal1;j<=$endsignal1;j++))
        do
         #./analysis /home/jf.fraga/events/ttbar_$j/Events/run_01/output_delphes.root /home/jf.fraga/analysis/temp/test_$j.root &
         ./analysis /disco1/SIMULACIONES/w+jets/w+jets_$j/Events/run_01/m_delphes_events.root /home/aflorez/SIMULATIONS/analysis/temp/test_$j.root &
	echo Running test_$j.root file
        done
        wait
	
####################################################################################

#Combine histograms

hadd -f /home/aflorez/SIMULATIONS/analysis/results/test.root /home/aflorez/SIMULATIONS/analysis/temp/test_*.root &
wait

#Removing temporal files
#cd /home/jf.fraga/temp/temporal/
#rm *
