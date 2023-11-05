set xlabel 'r'; set ylabel '{/Symbol r}'
plot [0:2] [-0.1:*] "results.dat" using 1:2 with lines title "{/Symbol r}(r)"
set terminal png
set output "density.png"