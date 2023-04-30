import numpy as np

mpiN = 0

step = 0.0005
x_ = np.arange(-0.003,0.003+step, step)
y_ = np.arange(-0.003,0.003+step, step)
mpi = [0,0.005, 0.01, 0.015, 0.025]
buffer = step * 3

z_ = np.arange(mpi[mpiN]-buffer,mpi[mpiN+1]+step, step)

print(mpi[mpiN])

u = 0
t = 884
ro = 22.65
p = 59.42e5
v = 1e-9

f = open("./configs/dime-1-4pr-cells/" + str(mpiN) +".dat", "w")

f.write("x	y	z	size	n[1]	n[2]	n[3]	vel	velX	velY	velZ	press	dens	temp\n")
sep = "\t"
for z in z_:
    for y in y_:
        for x in x_:
            f.write(str(x) + sep + str(y) + sep + str(z) + sep + str(v) + sep + "0"+ sep +"0"+ sep +"0"+ sep +"0"+ sep +"0"+ sep +"0"+ sep +"0"+ sep +str(p)+ sep +str(ro)+ sep +str(t) + "\n")

f.close()