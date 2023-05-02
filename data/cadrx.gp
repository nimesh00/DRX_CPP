# set autoscale xfix
# set autoscale yfix
# set autoscale cbfix
plot 'grid.dat' matrix nonuniform with image notitle
# show variables
bind "x" "exit gnuplot"
pause 0.3
reread