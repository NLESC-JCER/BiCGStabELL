 #Plot settings
set terminal pdf
set datafile separator ","
filename1="bicgstab_comparison_sleijpen.pdf"
filename2="bicgstab_comparison_tridiagonal.pdf"

#Make plot 1
set output filename1
set font 'Arial-Bold'
plot    'bicgstab_sleijpen.txt'  using 1:2 linestyle 1 with lines title 'BiCGStab',\
        'bicgstab2_sleijpen.txt'  using 1:2 linestyle 2 with lines title 'BiCGStab(2)'
set logscale y
set format y "10^{%T}"
set autoscale fix
set title 'Sleijpen example; -u_{xx}-u_{yy}+1000(xu_x+yu_y)+10u=f'
set xlabel 'Iteration' font 'Arial Bold,15'
set ylabel 'Residual ||Ax-b||_2' font 'Arial Bold,15'
set xtics font 'Arial Bold,15'
set ytics font 'Arial Bold,15'
set style line 1 lc rgb 'red' pt 5  pointsize 0.5 # square
set style line 2 lc rgb 'blue' pt 7 pointsize 0.5  # circle
set border linewidth 2
set key box linewidth 2 font 'Arial Bold,12'
set output filename1
replot

#Make plot 2
set output filename2
set font 'Arial-Bold'
plot    'bicgstab_tridiagonal_a_10.txt'  using 1:2 linestyle 1 with lines title 'BiCGStab',\
        'bicgstab2_tridiagonal_a_10.txt'  using 1:2 linestyle 2 with lines title 'BiCGStab(2)'
set logscale y
set format y "10^{%T}"
set autoscale fix
set title 'Tridiagonal example with a=10'
set xlabel 'Iteration' font 'Arial Bold,15'
set ylabel 'Residual ||Ax-b||_2' font 'Arial Bold,15'
set xtics font 'Arial Bold,15'
set ytics font 'Arial Bold,15'
set style line 1 lc rgb 'red' pt 5  pointsize 0.5 # square
set style line 2 lc rgb 'blue' pt 7 pointsize 0.5  # circle
set border linewidth 2
set key box linewidth 2 font 'Arial Bold,12'
set output filename2
replot
