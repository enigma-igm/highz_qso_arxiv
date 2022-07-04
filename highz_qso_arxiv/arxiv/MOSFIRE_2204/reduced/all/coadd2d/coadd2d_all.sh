for file in $ls *coadd2d.cfg;
do
    pypeit_coadd_2dspec --file $file
done