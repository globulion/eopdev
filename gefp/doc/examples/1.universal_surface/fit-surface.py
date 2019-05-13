#!/usr/bin/python3
"""
 Fit the universal surface of DMFT functional and plot it.
 Usage: ./fit-surface.py
 Note : modify input data within the script 
"""
from surface import fit_surface, plot_surface
def func_t_non(Z): return Z
def func_t_inv(Z): return 1./Z
def gen_gnuplot(data, surface_file):
    out = open('s.gpl', 'w')
    log = '#!/usr/bin/gnuplot\n'
    log+= 'set term eps\n'
    log+= "set output '%s'\n" % (surface_file[:-3] + 'eps')
    log+= "stats '%s' u (1./$2)\n" % surface_file
    for ia in data.pade.a:
        log += 'a%d%d = %14.7E\n' % tuple(ia)
    for ib in data.pade.b:
        log += 'b%d%d = %14.7E\n' % tuple(ib)
    log+= 'q(x,y) = ('
    loga= ''
    for ia in data.pade.a:
        loga += 'a%d%d * x**%f * y**%f +' % (ia[0], ia[1], ia[0], ia[1])
    log += loga[:-1] + ') / (1.0 + '
    logb = ''
    for ib in data.pade.b:
        logb+= 'b%d%d * x**%f * y**%f +' % (ib[0], ib[1], ib[0], ib[1])
    log += logb[:-1] + ')\n'
    log += 'set zrange[STATS_min:STATS_max]\n' 
    log += 'set view 30,300\n'
    log += "set xlabel 'Dynamic Correlation'\n"
    log += "set ylabel 'Non-Dynamic Correlation'\n"
    log += "set zlabel '1/t'\n"
    log += 'set ztics 20\n'
    log += "splot '%s' u (log($10/$1 + 1.0)):(log(2.0*$11/$1)):(1./$2) w p pt 7 ps 0.6 notit, q(x, y) notit\n"  % surface_file

    out.write(log)
    out.close()


# Input Data: Surface File and Function of Z
surface_file = 'test.dat'
func_t       = func_t_inv

# Input Data: 2D Pade Approximant Form
a_ord = [(0, 0), (0, 1), (1, 0), (1, 1)]
b_ord = [(0, 1), (1, 0), (1, 1)]

# Fit and Plot
data = fit_surface(surface_file, a_ord, b_ord, func_t=func_t, method='bfgs')
plot_surface(data)
print(" Fitted parameters:")
print(data.par)

t     = func_t(data.Z)
t_opt = func_t(data.Z_opt)

for i in range(len(t)):
    print("%14.5f %14.5f" % (t[i], t_opt[i]))

gen_gnuplot(data, surface_file)
