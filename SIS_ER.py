import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np
import os
import time

from networkx.utils import powerlaw_sequence

from ER import ER

class SISSimulation:
    prob_spontaneous_recovery = None
    prob_infection_range = [0]
    n_reps = None
    prob_being_initially_infected = None
    t_max = None
    t_trans = None
    output_path = None
    filename = None

    def __init__(self,
                 prob_spontaneous_recovery,
                 prob_infection_range,
                 n_reps,
                 prob_being_initially_infected,
                 t_max,
                 t_trans,
                 output_path=None):
        self.prob_spontaneous_recovery = prob_spontaneous_recovery

        for i in range(1, prob_infection_range):
            self.prob_infection_range.append(self.prob_infection_range[-1] + (1.0 / prob_infection_range))

        self.n_reps = n_reps
        self.prob_being_initially_infected = prob_being_initially_infected
        self.t_max = t_max
        self.t_trans = t_trans
        self.output_path = output_path

    def load_network(self, pajek_filepath, show=False):
        graph = nx.read_pajek(pajek_filepath)

        if show:
            nx.draw_circular(graph)
            plt.show()

        return graph

    def start(self):
        start_time = time.time()
        print '---starting simulation---'

        avrgs_of_infections = self.simulate_for_prob_infection_range(self.n_reps,
                                                                     self.prob_being_initially_infected,
                                                                     self.t_max,
                                                                     self.t_trans)

        self.plot(self.prob_infection_range, avrgs_of_infections)

        end_time = time.time() - start_time
        print '---simulation finished in {} hours, {} mins , {} secs ---'.format(round(end_time / 3600) % 24,
                                                                                 round(end_time / 60 % 60),
                                                                                 round(end_time % 60, 2))

    def simulate_for_prob_infection_range(self, n_reps, prob_being_initially_infected, t_max, t_trans):
        avrgs_of_infections = []
        network = ER().ER(500, 0.4)
        #network = nx.configuration_model(nx.utils.create_degree_sequence(n=500, exponent=2.7, sfunction=powerlaw_sequence))

        self.filename = os.path.join(self.output_path, 'A4' + 'ER_500_04_u01')
        nx.write_pajek(network, self.filename + '.net')

        adjacency_matrix = nx.to_numpy_matrix(network)
        infected_nodes = self.infecting_the_network(network, prob_being_initially_infected)
        for index, prob_infection in enumerate(self.prob_infection_range):
            avrgs_of_infections.append(self.simulate_for_prob_infection(network, infected_nodes, n_reps, prob_infection, t_max, t_trans, adjacency_matrix))
            print '\t{} infection probability ({}) = {}'.format(index, prob_infection, avrgs_of_infections[-1])

        return avrgs_of_infections

    # def simulate_for_prob_infection_range(self, n_reps, prob_being_initially_infected, t_max, t_trans):
    #     avrgs_of_infections = []
    #
    #     network = self.load_network('networks/airports_UW.net')
    #
    #     infected_nodes = self.infecting_the_network(network, prob_being_initially_infected)
    #     for index, prob_infection in enumerate(self.prob_infection_range):
    #         avrgs_of_infections.append(
    #             self.simulate_for_prob_infection(network, infected_nodes, n_reps, prob_being_initially_infected,
    #                                              prob_infection, t_max, t_trans))
    #         print '\t{} infection probability ({}) = {}'.format(index, prob_infection, avrgs_of_infections[-1])
    #
    #     return avrgs_of_infections

    def simulate_for_prob_infection(self, network, infected_nodes, n_reps, prob_infection, t_max, t_trans, adjacency_matrix):
        stationary = False
        avrg_over_repetitions = []
        for rep in range(1, n_reps + 1):
            avrg_over_time_steps = []
            for i in range(t_max, 0, -1):

                infected_nodes = self.get_next_state_partial_search(network, infected_nodes, prob_infection, adjacency_matrix)

                if stationary or i == (t_max - t_trans):
                    avrg_over_time_steps.append(len(infected_nodes) / float(len(network.nodes())))
                    stationary = True

            avrg_over_repetitions.append(sum(avrg_over_time_steps) / float(len(avrg_over_time_steps)))

        return sum(avrg_over_repetitions) / float(len(avrg_over_repetitions))

    def get_next_state_partial_search(self, network, infected_nodes, prob_infection, adjacency_matrix):
        infected_nodes_next = set()
        susceptible_nodes_next = set()
        visited_nodes = np.zeros(len(network.nodes()))
        for node in infected_nodes:
            if random.uniform(0, 1) < self.prob_spontaneous_recovery:
                susceptible_nodes_next.add(node)

            neighbors = set(network.neighbors(node)).difference(infected_nodes)

            for susc_node in neighbors:
                if visited_nodes[susc_node] == 1:
                    continue

                visited_nodes[susc_node] = 1
                infected_found = 0

                susc_neighbors = network.neighbors(susc_node)
                if len(susc_neighbors) > len(infected_nodes):
                    for inf in infected_nodes:
                        if adjacency_matrix.item((susc_node, inf)) == 1:
                            infected_found += 1
                else:
                    for susc_neighbor in susc_neighbors:
                        if susc_neighbor in infected_nodes:
                            infected_found += 1

                while infected_found > 0:
                    if random.uniform(0, 1) < prob_infection:
                        infected_nodes_next.add(susc_node)
                        break
                    infected_found -= 1

        return infected_nodes.difference(susceptible_nodes_next).union(infected_nodes_next)

    def infecting_the_network(self, graph, p_being_initially_infected):
        infected_nodes = set()

        while True:
            for node in graph.nodes():
                if random.uniform(0, 1) < p_being_initially_infected:
                    infected_nodes.add(node)

            if len(infected_nodes) > 0:
                break

        return infected_nodes

    def create_pajek_file(self, network):
        network = nx.write_pajek(network, self.output_path)

    def plot(self, x, y):
        title = 'ER (500, 0.4)'
        title += '\nSIS (' + r'$\mu = {},\rho$0 = {})'.format(self.prob_spontaneous_recovery, self.prob_being_initially_infected)
        title += '\nN-rep={},   T-max={},   T-trans={}'.format(self.n_reps, self.t_max, self.t_trans)
        plt.title(title)
        plt.xlabel(r'$\beta$', fontsize=18)
        plt.ylabel(r'$\rho$', fontsize=18)

        plt.plot(x, y, linestyle='-')
        axes = plt.gca()
        axes.set_xlim([0, 1])
        axes.set_ylim([0, 1])
        plt.grid()
        plt.tight_layout()

        plt.savefig(self.filename + '.png')
        plt.clf()

if __name__ == '__main__':
    sis = SISSimulation(0.5, 51, 50, 0.2, 1000, 900, 'outputs/')
    sis.start()