ls output_folder > list_folder.dat
for x in `cat list_folder.dat` ; do
rm output_folder/$x/* -Rf
done
rm list_folder.dat
