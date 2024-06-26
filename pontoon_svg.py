# create a nacaxx svg file for pontoon
import svgwrite
from aerofoil import naca_foil

# make a naca foil
x,y = naca_foil(0.25,0)

# move the foil to the right
# x = [i+0.5 for i in x]
# # move the foil down
# y = [i+0.5 for i in y]

# scale to
size = 1800 # cm
x = [i*size for i in x]
y = [i*size for i in y]
z = [0 for i in x]


# make an svg file with the foil coordinates
dwg = svgwrite.Drawing('pontoon.svg', profile='tiny')
dwg.add(dwg.polyline(points=[(x[i],y[i]) for i in range(len(x))], fill='none', stroke='black'))
dwg.save()

# save the foil coordinates to a csv file
f = open('pontoon.csv','w')
for i in range(len(x)):
    f.write(str(x[i])+','+str(y[i])+','+str(z[i])+'\n')
f.close()

# save the foil coordinates to a txt file
f = open('pontoon.txt','w')
for i in range(len(x)):
    f.write(str(x[i])+' '+str(y[i])+' '+str(z[i])+'\n')
f.close()



