import os
import networkx as nx
import matplotlib.pyplot as plt
import random
from scipy import stats

import numpy as np

from A2.ER import ER
from plotter.Plotter import Plotter

class BA:
    def BA(self, previous_nodes, n, m, output_path=None):
        if (m < 1) or (m > previous_nodes):
            raise ValueError('\nFacts. Must. \nm >= 1 (' + str(m >= 1) + ') and\nprevious_nodes >= m (' + str(previous_nodes >= m) + ') \n')

        G = ER().ER(previous_nodes, 1)  # start with a network of "previous_nodes" nodes full connected

        begin_nodes = G.nodes()

        for node in range(0, n):
            total_edges = sum(nx.degree(G).values())
            degree_values = nx.degree(G)
            degree_segmented_values = {}
            next_node = G.number_of_nodes()

            for key, val in degree_values.iteritems():
                if key == 0:
                    degree_segmented_values[key] = [0.0, val - 0.1]
                else:
                    degree_segmented_values[key] = [(degree_segmented_values.values()[-1][1] + 0.1),
                                                    (degree_segmented_values.values()[-1][1] + val)]

            previous_selected_node = []
            count = 0
            while count < m:
                random_number = random.uniform(0.0, 1.0) * total_edges

                for key, val in degree_segmented_values.iteritems():
                    if (random_number >= val[0]) and (random_number <= val[1]):
                        if key not in previous_selected_node:
                            G.add_edge(key, next_node)
                            #print 'edge (' + str(key) + ', ' + str(next_node) + ')'
                            previous_selected_node.append(key)
                            count += 1
                        break

        end_nodes = set(G.nodes()) - set(begin_nodes)

        if output_path is not None:
            title = 'Barabasi & Albert G = (n, m)  -->  G = ({}, {})  --> Previous nodes = {}'.format(n, m, previous_nodes)
            title = title + '\nEstimation of the exponent is = ' + str(self.getEstimationOfTheExponent(G))
            plt.title(title)
            nx.draw_networkx_labels(G, pos=nx.shell_layout(G), font_size=10, font_family='sans-serif')
            nx.draw_networkx_nodes(G, pos=nx.shell_layout(G), nodelist=begin_nodes, node_color='r')
            nx.draw_networkx_nodes(G, pos=nx.shell_layout(G), nodelist=end_nodes, node_color='g')
            nx.draw_networkx_edges(G, pos=nx.shell_layout(G))
            plt.savefig(os.path.join(output_path, 'BA' + '_' + str(previous_nodes) + '_' + str(n) + '_' + str(m) + '.png'))
            plt.clf()

            adjacency_matrix = nx.to_numpy_matrix(G)
            max_k = previous_nodes + n
            binomial_distribution = self.getProbDegreeDistribution(max_k, m)

            plotter = Plotter()
            plotter.histLogLogScalePDF_BA(binomial_distribution, previous_nodes, n, m, adjacency_matrix, output_path)

        return G

    def getProbDegreeDistribution(self, max_k, m):
        dist = []

        for i in range(0, m):
            dist.append(0.0)

        for k in range(m, max_k):
            dist.append(((2 * m) * (m + 1)) / (k * (k + 1) * (k + 2.0)))

        return dist

    def getEstimationOfTheExponent(self, graph):
        degree_log = [np.math.log10(x) for x in graph.degree().values()]

        k_min_log = min(degree_log)
        k_max_log = max(degree_log) + 1

        # calculating the range between bins
        bins = 10
        binsMargin = [k_min_log]
        binsRange = (k_max_log - k_min_log) / bins
        for i in range(0, bins):
            binsMargin.append(binsMargin[i] + binsRange)

        #
        intervals = [0] * len(binsMargin)
        for i in range(0, len(binsMargin)):
            for k_log in degree_log:
                if k_log >= binsMargin[i] and k_log < binsMargin[i + 1]:
                    intervals[i] = intervals[i] + 1

        # getting probabilities for each bin
        bins_prob = [x / float(len(degree_log)) for x in intervals]

        # linear regression
        gradient, intercept, r_value, p_value, std_err = stats.linregress(binsMargin, bins_prob)

        return -gradient

    def getArrayValueRepetitions(self, data):
        y = np.bincount(data)
        ii = np.nonzero(y)[0]
        return dict(zip(ii, y[ii]))

if __name__ == '__main__':
    generator = BA()

    #BA: different values of "m" (number of edges that each new nodes forms with the existing nodes), e.g. m=1, 2, 4, 10
    generator.BA(10, 10, 1, 'YOUR_OUTPUT_HERE')