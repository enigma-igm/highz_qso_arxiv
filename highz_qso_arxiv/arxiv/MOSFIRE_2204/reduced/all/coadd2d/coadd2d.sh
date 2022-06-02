for z in `ls *.cfg`
do
    # echo $z
    pypeit_coadd_2dspec --file $z
done
