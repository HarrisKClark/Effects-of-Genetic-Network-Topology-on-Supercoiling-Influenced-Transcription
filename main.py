import random
from Polymerase import Polymerase
from genome_handler import Genome
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

genome_size = 1600 #bp

gene1_size = 100 #bp
gene1_position = 400
gene1_direction = -1 #clockwise
gene1_time = 0

gene2_size = 100 #bp
gene2_position = 800
gene2_direction = 1 #clockwise
gene2_time = 400

gene3_size = 100 #bp
gene3_position = 1200 #bp
gene3_direction = -1 #clockwise
gene3_time = 800

gene4_size = 100 #bp
gene4_position = 1600 #bp
gene4_direction = -1 #clockwise
gene4_time = 0

gene5_size = 100 #bp
gene5_position = 2000 #bp
gene5_direction = 1 #clockwise
gene5_time = 0

gene6_size = 100 #bp
gene6_position = 2400 #bp
gene6_direction = -1 #clockwise
gene6_time = 0

gene7_size = 100 #bp
gene7_position = 2800 #bp
gene7_direction = 1 #clockwise
gene7_time = 0

gene8_size = 100 #bp
gene8_position = 3200 #bp
gene8_direction = -1 #clockwise
gene8_time = 0

gene1 = [gene1_size, gene1_position, gene1_direction, gene1_time]
gene2 = [gene2_size, gene2_position, gene2_direction, gene2_time]
gene3 = [gene3_size, gene3_position, gene3_direction, gene3_time]
gene4 = [gene4_size, gene4_position, gene4_direction, gene4_time]
gene5 = [gene5_size, gene5_position, gene5_direction, gene5_time]
gene6 = [gene6_size, gene6_position, gene6_direction, gene6_time]
gene7 = [gene7_size, gene7_position, gene7_direction, gene7_time]
gene8 = [gene8_size, gene8_position, gene8_direction, gene8_time]

genes=[gene1,gene2, gene3]
#genes = [gene1, gene2, gene3, gene4]
#genes = [gene1,gene2, gene3, gene4,gene5,gene6, gene7,gene8]

genome = Genome(genome_size, genes)

def global_SC(polys,z):
    '''
    returns the global SC function at a paticular z, given an array of polymerases

    0.06 for base level SC
    '''

    SC = 0

    for poly in polys:
        SC += poly.local_SC(z)


    return (SC+ 0.06)


def topo_global_SC(polys, z):
    '''
    returns the global SC function at a paticular z, given an array of polymerases

    0.06 for base level SC
    '''

    a = 0.8

    SC = 0

    for poly in polys:
        SC += poly.local_SC(z)

    return (SC + 0.06)*a

'''
simulation rules:

check to add polymerases

move polymerases

if polyz > coding length remove poly, then add mrna

'''

#first pass

polymerases = []
global_SC0 = 0.06  # base level, to use when no polys attatched
promoter_SC = [0, 0, 0]
mRNA = [[0] for _ in range(len(genes))]
delm = 0.1

tspan = 4000

yvals= []

for i in range(1,tspan):
    percent = i/10
    print(percent)

    yvalnew = []
    for j in range(1,genome_size+1):
        yvalnew.append(global_SC(polymerases, j))

    yvals.append(yvalnew)

    if polymerases == []:
        promoter_SC = [global_SC0]*len(genes)
    else:
        promoter_SC=[]
        for j in range(0,len(genes)):
            promoter_SC.append(global_SC(polymerases, genome.genes[j][1]))

    add_poly = genome.add_polymerase_new(promoter_SC, i)

    if add_poly != -1:
        new_poly = Polymerase(genome.genes[add_poly][2], genome.genes[add_poly][1], add_poly) # direction and position
        polymerases.append(new_poly)

    for poly in polymerases:
        poly.z += poly.dzdt(global_SC(polymerases, poly.z))

    for i in range(0,len(mRNA)):
        mRNA[i].append(mRNA[i][-1])

    #remove poly + add mRNA
    for poly in polymerases:
        if abs(genome.genes[poly.gene][1] - poly.z) > genome.genes[poly.gene][0]: #gone through the whole thing
            polymerases.pop(polymerases.index(poly)) #remove polymerase

            mRNA[poly.gene][-1]+=1

    #remove mRNA

    for i in range(0,len(mRNA)):
        remove_prop = 0.001*mRNA[i][-1]
        add_prop = 0.1

        tot = abs(remove_prop+add_prop)

        propnum = random.uniform(0,1)*tot

        if propnum < remove_prop:
            mRNA[i][-1] -=1


t = np.linspace(1,tspan,tspan)


for i in range(0,len(genes)):
    label = "gene: "+ str(i)
    if i ==0:
        plt.plot(t, mRNA[i], label=label, linewidth=2,color="blue")
    elif i==1:
        plt.plot(t, mRNA[i], label=label, linewidth=2, color="red")
    elif i==  2:
        plt.plot(t, mRNA[i], label=label, linewidth=2, color="green")
    elif i==  3:
        plt.plot(t, mRNA[i], label=label, linewidth=2, color="gold")
    elif i==  4:
        plt.plot(t, mRNA[i], label=label, linewidth=2, color="orange")
    elif i==  5:
        plt.plot(t, mRNA[i], label=label, linewidth=2, color="purple")
    elif i==  6:
        plt.plot(t, mRNA[i], label=label, linewidth=2, color="gray")
    else:
        plt.plot(t, mRNA[i], label=label, linewidth=2, color="brown")

#plt.plot(t,mRNA[0], label="gene 1")
#plt.plot(t,mRNA[1], label="gene 2")
#plt.plot(t,mRNA[2], label="gene 3")
#plt.plot(t,mRNA[3], label="gene 4")
plt.legend()
plt.show()


# Create an animation of yvals at each time t
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)

ax.set_xlim(1, genome_size)
# Determine y-limits based on yvals data
y_min = min(min(yval) for yval in yvals)
y_max = max(max(yval) for yval in yvals)
ax.set_ylim(y_min, y_max)
ax.set_xlabel('Genome Position')
ax.set_ylabel('Global SC')
ax.set_title('Global Supercoiling Over Time')

# Plot gene positions as vertical lines
for gene in genes:
    ax.axvline(x=gene[1], color='gray', linestyle='--')

# Initialization function
def init():
    line.set_data([], [])
    return line,

# Animation function
def animate(i):
    x = np.linspace(1, genome_size, genome_size)
    y = yvals[i]
    line.set_data(x, y)
    return line,

# Create the animation
ani = animation.FuncAnimation(
    fig, animate, init_func=init,
    frames=len(yvals), interval=50, blit=True
)

#writervideo = animation.PillowWriter(fps=120)
#ani.save('filename.gif', writer=writervideo)
#plt.close()
print("done")

'''
plt.plot(np.linspace(1,1800,1800),yval)
plt.plot(np.linspace(1,1800,1800), [0]*1800, "k")
plt.axvline(x=gene1_position)
plt.axvline(x=gene2_position)
plt.axvline(x=gene3_position)
plt.show()
'''
