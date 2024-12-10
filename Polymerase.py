import math

class Polymerase:
    def __init__(self, direction, z, gene):
        self.optimal_SC = -0.2  # optimal SC density

        self.direction = direction #direction of travel
        self.z = z #current position of polymerase
        self.gene = gene

    def dzdt(self, global_SC):
        '''
        returns the velocity of the polymerase (bp/s)

        uses global_SC @ self.z - to tune speed

        optimal speed (@ global_sc) = 60bp/s
        '''

        optimal_speed = 6

        return optimal_speed * self.direction*math.e**(-(global_SC-self.optimal_SC)**2)/0.5

    def local_SC(self, z):
        '''
        returns valye of the local supercoiling function @ position z
        '''

        a = 0.1 # supercoiling addition constant
        b = 100 # large number for sigmoid function

        try:
            return a * self.direction * (1 / (1 + math.e ** (-b * (z - self.z))) - 0.5)
        except OverflowError:
            return a * self.direction * (- 0.5)
