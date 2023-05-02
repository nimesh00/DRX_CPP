set autoscale;

filename1 = "../data/0.01.csv"
filename2 = "../data/strainVstress.dat"

plot filename1 using 1:2 with lines, filename2 using 1:2 with lines

# bind "x" "exit gnuplot"
