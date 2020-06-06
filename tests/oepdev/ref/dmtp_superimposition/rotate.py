#!/usr/bin/python2.7
import dma

dma_1 = dma.DMA('w1.par')
dma_2 = dma.DMA('w2.par')

# task: rotate dma_1 onto dma_2
rms = dma_1.sup(dma_2.get_pos())
dma_1.MAKE_FULL()
p = dma_1.DMA_FULL
print(dma_1)

c, q, d, Q, O, H = p

label = {0:"X", 1:"Y", 2:"Z"}

def print_dipole(d):
    log = "const double d_ref[%d] =/\n{" % d.size
    for i in range(len(d)):
        di = d[i]
        log += "%13.6E,"*3 % tuple(di) + "\n"
    log = log[:-2] + "};"
    print(log)

def print_quadrupole(d):
    log = "const double Q_ref[%d] =/\n{" % d.size
    for i in range(len(d)):
        di = d[i]
        for a in range(3):
            for b in range(3):
              if b>=a:
                l = "Q%c%c" % (label[a], label[b])
                #print(l)
                log += "%13.6E," % di[a,b]
    log = log[:-1] + "};"
    print(log)

def print_octupole(d):
    log = "const double O_ref[%d] =/\n{" % d.size
    for i in range(len(d)):
        di = d[i]
        for a in range(3):
            for b in range(3):
              if b>=a:
               for c in range(3):
                if c>=b:
                 l = "O%c%c%c" % (label[a], label[b], label[c])
                 #print(l)
                 log += "%13.6E," % di[a,b,c]
    log = log[:-1] + "};"
    print(log)

def print_hexadecapole(d):
    log = "const double H_ref[%d] =/\n{" % d.size
    for i in range(len(d)):
        di = d[i]
        for a in range(3):
            for b in range(3):
              if b>=a:
               for c in range(3):
                if c>=b:
                 for u in range(3):
                  if u>=c:
                    l = "H%c%c%c%c" % (label[a], label[b], label[c], label[u]) 
                    #print(l)
                    log += "%13.6E," % di[a,b,c,u]
    log = log[:-1] + "};"
    print(log)

def symmetry_labels(*labels):
    labels = list(labels)
    labels.sort()
    return labels
    
def make_code_for_hexadecapoles():
    log = ""
    for a in range(3):
      for b in range(3):
       if b>=a:
         for c in range(3):
          if c>=b:
            for d in range(3):
             if d>=c:
              new_h_label = "rH" + label[a] + label[b] + label[c] + label[d]
              mini_log = ""
              for i in range(3):
               ria = "R" + label[i] + label[a]
               for j in range(3):
                rjb = "R" + label[j] + label[b]
                for k in range(3): 
                 rkc = "R" + label[k] + label[c]
                 for l in range(3):
                  rld = "R" + label[l] + label[d]
                  ii, jj, kk, ll = symmetry_labels(i,j,k,l)
                  h_label = "H" + label[ii] + label[jj] + label[kk] + label[ll]
                  mini_log += "                %s * %s * %s * %s * %s +\n" % (ria, rjb, rkc, rld, h_label)
              mini_log = mini_log[:-2] + ";\n"
              log += "double %s = \n%s" % (new_h_label, mini_log)
    print log

print_dipole(c)
print_dipole(d)
print_quadrupole(Q)
print_octupole(O)
print_hexadecapole(H)

make_code_for_hexadecapoles()    


#print(dma_1-dma_2)
print(rms)
