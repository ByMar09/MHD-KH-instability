reset
set terminal pngcairo  enhanced color size 1500, 2300	
set output '00_1d_test.png' 
set encoding utf8
set view map
set border 4095 front linetype -1 linewidth 1.000
#unset surface
unset key


######## README ###########

# To use this script you need a merge file with data of the two simulations that you would to compare
# to do this you could use the comand
# paste -d" " file1.dat file2.dat > file3.dat | sort
# to genearte a file3.dat with data from print_backup_xxxx.dat type files that contains the adequate format to print


#margins takes four numbers set multiplot margins <left>,<right>,<bottom>,<top>,
#which give the fixed overall margins around the multiplot layout.
#spacing takes two number set multiplot spacing <xspacing>,<yspacing>
#which give the distance between two rows (<yspacing>) or two columns (<xspacing>).


######## README ###########

### Start multiplot (3x3 layout)

set multiplot layout 5,3 

#set lmargin 6
#set rmargin 100006
#set tmargin 1
#set bmargin 1

unset xlabel
#set xrange [ -6.0 : 6.0 ] noreverse nowriteback
#set xtics    -6.0, 3.0, 6.0
#set yrange [ -6.0 : 6.0 ] noreverse nowriteback
#set ytics    -6.0, 3.0, 6.0
#set zlabel "P" font "arial,10"
#set zrange [ -6.000 : 6.000  ] noreverse nowriteback 
#set cblabel ""  font "arial,10"
set colorbox vertical

set format cb "%.1e"


# --- GRAPH (1,1)

#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set format x ""

set label 1 'Bx' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($25) w l


# --- GRAPH (1,2)


#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

unset xlabel
unset ylabel

set label 1 'By ' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($26) w l


# --- GRAPH (1,3)


#set cbrange [  0 : 1  ]
#set cbtics     0, 0.5 , 1

unset title
unset xlabel
unset ylabel 


set label 1 'Bz' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($27) w l


# --- GRAPH (2,1)


unset xlabel
set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced


#set cbrange [  -2.5e-5 : 2.5e-5  ]
#set cbtics     -2.e-5, 1e-5 , 2.e-5

#set format cb "%.1t{/Symbol \264}10^{%L}" 
#set cbtics add ("0" 0) 

set label 1 'Vx' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($28) w l

# --- GRAPH (2,2)



#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

unset title
unset ylabel 



set label 1 'Vy' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($29) w l

# --- GRAPH (2,3)


#set cbrange [  -1 : 1  ]
#set cbtics     -1, 0.5 , 1

unset ylabel 

set label 1 'Vz' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($30) w l

# --- GRAPH (3,1)


#set cbrange [  0.49 : 0.57  ]
#set cbtics     0.49, 0.02 , 0.57

set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set label 1 'Ex' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($31) w l

# --- GRAPH (3,2)


#set cbrange [  1 : 1.014  ]
#set cbtics     1, 0.002 , 1.014

unset ylabel 

set label 1 'Ey' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($32) w l


# --- GRAPH (3,3)


#set cbrange [  0.5 : 0.514  ]
#set cbtics     0.5, 0.002 , 0.514

unset ylabel 

set label 1 'Ez' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($33) w l

# --- GRAPH (4,1)

#set cbrange [  0 : 120  ]
#set cbtics     0, 30 , 120

set ylabel "y" font "Times-Roman,20"  rotate by 0 enhanced

set label 1 'Q' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($34) w l

# --- GRAPH (4,2)

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

unset ylabel 

set label 1 'Rho ' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($35) w l

# --- GRAPH (4,3)

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

set format x "%g"

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 

set label 1 'P ' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($36) w l

# --- GRAPH (5,1)

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 

set label 1 'psi-->(divE) ' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($23) w l

# --- GRAPH (5,2)

#set cbrange [  0 : 4.5e-4  ]
#set cbtics     0, 2e-4 , 4e-4

set xlabel "x" font "Times-Roman,20"  rotate by 0 enhanced
unset ylabel 

set label 1 'phi-->(divB) ' at graph 0.05,0.85 font ',14' front
plot "file3.dat" u ($22):($24) w l


unset multiplot
### End multiplot

