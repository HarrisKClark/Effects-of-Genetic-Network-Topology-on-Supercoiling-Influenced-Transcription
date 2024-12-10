import math
import random

class Genome:
    def __init__(self, genome_size, genes):
        self.optimal_SC = -0.2  # optimal SC density

        self.genome_size = genome_size
        self.genes = genes

    def poly_add_energy(self, SC):
        '''
        returns the binding energy given a supercoiling value [0-1]

        1 at optimal SC

        s is tunable to KbT and describes binding sensitivity to SC
        '''

        s = 0.2

        return math.e**(-((SC-self.optimal_SC)**2)/s)

    def add_polymerase(self, SCs):
        '''
        decides if and where to add a polymerase based of the global supercoiling
        function values at the promoter location

        SCs is an array of the SC value at each promoter
        '''

        binding_position = -1

        max_energy = len(self.genes) #max valye of the enrgy func is 1

        binding_energy1 = self.poly_add_energy(SCs[0])
        binding_energy2 = self.poly_add_energy(SCs[1])
        binding_energy3 = self.poly_add_energy(SCs[2])

        random_energy = max_energy*random.uniform(0, 1)

        if random_energy <= binding_energy1:
            binding_position = 0
        elif random_energy <= (binding_energy1 + binding_energy2):
            binding_position = 1
        elif random_energy <= (binding_energy1 + binding_energy2 + binding_energy3):
            binding_position = 2

        return binding_position

    def add_polymerase_new(self, SCs, t):
        '''
        decides if and where to add a polymerase based of the global supercoiling
        function values at the promoter location

        SCs is an array of the SC value at each promoter
        '''
        binding_position = -1
        max_energy = len(self.genes)
        binding_energies = []
        null = 1

        for i in range(0,len(SCs)):
            if self.genes[i][3] > t:
                null = 0
            binding_energies.append(null*self.poly_add_energy(SCs[i]))
            null = 1

        random_energy = max_energy * random.uniform(0, 1)

        total_engr = 0

        for i in range(0, len(SCs)):
            total_engr+=binding_energies[i]

            if random_energy <= total_engr:
                binding_position = i
                break

        return binding_position
