import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np
import os
import time

from A2.BA import BA
from A2.ER import ER


class SISSimulation:
    prob_spontaneous_recovery = 0.0
    prob_infection_range = [0]
    output_path = None

    def __init__(self, prob_spontaneous_recovery, prob_infection_range, output_path=None):
        self.prob_spontaneous_recovery = prob_spontaneous_recovery

        for i in range(1, prob_infection_range):
            self.prob_infection_range.append(self.prob_infection_range[-1] + (1.0 / prob_infection_range))

        self.output_path = output_path

    def load_network(self, pajek_filepath, show=False):
        graph = nx.read_pajek(pajek_filepath)

        if show:
            nx.draw_circular(graph)
            plt.show()

        return graph

    def start(self, n_reps, prob_being_initially_infected, t_max, t_trans):
        start_time = time.time()
        print '---starting simulation---'

        avrgs_of_infections = self.simulate_for_prob_infection_range(n_reps, prob_being_initially_infected, t_max, t_trans)

        self.plot(self.prob_infection_range, avrgs_of_infections)

        end_time = time.time() - start_time
        print '---simulation finished in {} hours, {} mins , {} secs ---'.format(3600 % 24, round(end_time / 60 % 60), round(end_time % 60, 2))

    def simulate_for_prob_infection_range(self, n_reps, prob_being_initially_infected, t_max, t_trans):
        avrgs_of_infections = []

        for index, prob_infection in enumerate(self.prob_infection_range):
            avrgs_of_infections.append(self.simulate_for_prob_infection(n_reps, prob_being_initially_infected, prob_infection, t_max, t_trans))
            print '\t{} infection probability ({}) = {}'.format(index, prob_infection, avrgs_of_infections[-1])

        return avrgs_of_infections

    def simulate_for_prob_infection(self, n_reps, prob_being_initially_infected, prob_infection, t_max, t_trans):
        network = ER().ER(50, 0.1)
        #network = BA().BA(500, 10, 3)

        infected_nodes = self.infecting_the_network(network, prob_being_initially_infected)
        #print 'initial infected nodes => {} / {} = {}'.format(len(infected_nodes), float(len(network.nodes())), len(infected_nodes) / float(len(network.nodes())))

        stationary = False

        avrg_over_repetitions = []

        for rep in range(1, n_reps + 1):
            #print '\nRep #{}'.format(rep)
            avrg_over_time_steps = []

            for i in range(t_max, 0, -1):
                infected_nodes_next = set()
                susceptible_nodes_next = set()

                for node in network.nodes():
                    if node in infected_nodes:
                        if random.uniform(0, 1) < self.prob_spontaneous_recovery:
                            susceptible_nodes_next.add(node)
                    else:
                        neighbors = network.neighbors(node)

                        if len(neighbors) > 0:
                            infected_neighbors = 0

                            for neighbor in neighbors:
                                if neighbor in infected_nodes:
                                    infected_neighbors += 1

                            while infected_neighbors > 0:
                                if random.uniform(0, 1) < prob_infection:
                                    infected_nodes_next.add(node)
                                    break
                                infected_neighbors -= 1


                infected_nodes = infected_nodes.difference(susceptible_nodes_next).union(infected_nodes_next)
                #print len(infected_nodes)

                if stationary or i == (t_max - t_trans):
                    #print 'step #{} avrg_over_time -> {} / {} = {}'.format(i, len(infected_nodes), float(len(network.nodes())), len(infected_nodes) / float(len(network.nodes())))
                    avrg_over_time_steps.append(len(infected_nodes) / float(len(network.nodes())))
                    stationary = True

            avrg_over_repetitions.append(sum(avrg_over_time_steps) / float(len(avrg_over_time_steps)))
            #print 'avrg_over_repetitions = {}'.format(sum(avrg_over_time_steps) / float(len(avrg_over_time_steps)))

        return sum(avrg_over_repetitions) / float(len(avrg_over_repetitions))

    def infecting_the_network(self, graph, p_being_initially_infected):
        infected_nodes = set()

        while True:
            for node in graph.nodes():
                if random.uniform(0, 1) < p_being_initially_infected:
                    infected_nodes.add(node)

            if len(infected_nodes) > 0:
                break

        return infected_nodes

    def plot(self, x, y):
        plt.plot(x, y, linestyle='-')
        axes = plt.gca()
        axes.set_xlim([0, 1])
        axes.set_ylim([0, 1])
        plt.grid()
        plt.savefig(os.path.join(self.output_path, 'A4' + '_' + str(random.randint(0, 10000)) + '.png'))
        plt.clf()

if __name__ == '__main__':
    sis = SISSimulation(1.0, 51, 'YOUR_OUTPUT_PATH')
    sis.start(10, 0.2, 1000, 900)