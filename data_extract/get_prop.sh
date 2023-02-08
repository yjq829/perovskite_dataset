#!/bin/bash
B=`pwd`
for ((  n = 1;  n <= 550;  n++  ))
do
A=`head -$n sys | tail -1 | cut -f 1`
C=`head -$n s1 | tail -1 | cut -f 1`
D=`head -$n s2 | tail -1 | cut -f 1`


##  Analyize PBE Relax Calc  ##

#cp $B/band_edges.py $B/$A/PBE_relax/
#cp $B/my_colours.conf $B/$A/PBE_relax/

#cd $B/$A/PBE_relax/
#sumo-dosplot --config my_colours.conf
#python band_edges.py
#cd $B

#sed -n '1p' $B/$A/PBE_relax/gap.txt > $B/$A/PBE_relax/gap
#sed -n '2p' $B/$A/PBE_relax/gap.txt > $B/$A/PBE_relax/dir_gap
#sed -n '3p' $B/$A/PBE_relax/gap.txt > $B/$A/PBE_relax/dos_gap
#cat $B/$A/PBE_relax/gap >> $B/gap
#cat $B/$A/PBE_relax/dir_gap >> $B/dir_gap
#cat $B/$A/PBE_relax/dos_gap >> $B/dos_gap

#cat $B/$A/PBE_relax/vbm_bs.txt $B/x >> $B/vbm
#cat $B/$A/PBE_relax/cbm_bs.txt $B/x >> $B/cbm
#cat $B/$A/PBE_relax/gap.txt >> $B/gap
#cat $B/$A/PBE_relax/aa.txt >> $B/a
#cat $B/$A/PBE_relax/bb.txt >> $B/b
#cat $B/$A/PBE_relax/cc.txt >> $B/c
#cat $B/$A/PBE_relax/vol.txt >> $B/vol
#cat $B/$A/PBE_relax/alpha.txt >> $B/alpha
#cat $B/$A/PBE_relax/beta.txt >> $B/beta
#cat $B/$A/PBE_relax/gamma.txt >> $B/gamma
#cat $B/$A/PBE_relax/toten.txt >> $B/toten

#rm $B/$A/PBE_relax/gap_shell
#cd $B/$A/PBE_relax/
#$B/vtstscripts-935/bandgap.pl > gap_info
#grep "direct band gap" gap_info >> gap_shell
#cd $B
#cat $B/$A/PBE_relax/gap_shell >> $B/gaps_shell


##  Analyze LOPTICS Calc  ##

#rm $B/$A/LOPTICS/REAL.in
#rm $B/$A/LOPTICS/eps*
#rm $B/$A/LOPTICS/fom*
#cp $B/Data_AM.xlsx $B/$A/LOPTICS/
#cp $B/get_fom.py $B/$A/LOPTICS/
#cp $B/optics.sh $B/$A/LOPTICS/

#cd $B/$A/LOPTICS/
#sumo-optplot --anisotropic --height 6 --width 8 --xmax 8
#python get_fom.py
#rm Data_AM.xlsx
#./optics.sh &
#cd $B

head -1 $B/$A/LOPTICS/REAL.in > $B/$A/LOPTICS/eps
cat $B/$A/LOPTICS/eps >> $B/eps
cat $B/$A/LOPTICS/fom.txt $B/x >> $B/fom


#cp $B/am1.5G.dat $B/$A/LOPTICS/
#cp $B/ss.py $B/$A/LOPTICS/
#echo $C > $B/$A/LOPTICS/s1
#echo $D > $B/$A/LOPTICS/s2
#cd $B/$A/LOPTICS/
#cat s1 s2 ss.py >> calc_slme.py
#rm s1 s2 ss.py
#python calc_slme.py
#cd $B

#sed -n '1p' $B/$A/LOPTICS/slme_data.txt > $B/$A/LOPTICS/s1
#sed -n '2p' $B/$A/LOPTICS/slme_data.txt > $B/$A/LOPTICS/s2
#sed -n '3p' $B/$A/LOPTICS/slme_data.txt > $B/$A/LOPTICS/s3
#cat $B/$A/LOPTICS/s1 >> $B/slme_5um
#cat $B/$A/LOPTICS/s2 >> $B/slme_100um
#cat $B/$A/LOPTICS/s3 >> $B/sat_thick



##  Analyize HSE Relax Calc  ##

#cp $B/band_edges.py $B/$A/HSE_relax_2x2x2/
#cp $B/my_colours.conf $B/$A/HSE_relax_2x2x2/

#cd $B/$A/HSE_relax_2x2x2/
#sumo-dosplot --config my_colours.conf
#python band_edges.py
#cd $B

#cat $B/$A/HSE_relax_2x2x2/vbm_bs.txt $B/x >> $B/vbm_hse
#cat $B/$A/HSE_relax_2x2x2/cbm_bs.txt $B/x >> $B/cbm_hse
#cat $B/$A/HSE_relax_2x2x2/gap.txt >> $B/gap_hse
#cat $B/$A/HSE_relax_2x2x2/aa.txt >> $B/a_hse
#cat $B/$A/HSE_relax_2x2x2/bb.txt >> $B/b_hse
#cat $B/$A/HSE_relax_2x2x2/cc.txt >> $B/c_hse
#cat $B/$A/HSE_relax_2x2x2/vol.txt >> $B/vol_hse
#cat $B/$A/HSE_relax_2x2x2/alpha.txt >> $B/alpha_hse
#cat $B/$A/HSE_relax_2x2x2/beta.txt >> $B/beta_hse
#cat $B/$A/HSE_relax_2x2x2/gamma.txt >> $B/gamma_hse
#cat $B/$A/HSE_relax_2x2x2/toten.txt >> $B/toten_hse



done
