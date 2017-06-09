import os
import networkx as nx
import matplotlib.pyplot as plt
import random
import math

from plotter.Plotter import Plotter

class ER:
    """Erdos-Renyi G = (n, p).

    ER: different values of "K" for G(N,K), or of "p" for G(N,p), e.g. p=0.00, 0.01, 0.02, 0.03, 0.05, 0.1
    Primero se crea un grafo de "n" nodos, luego se van tomando parejas de nodos aleatorias que no tengan aun enlace
    entre ellos y se toma un valor numero aleatorio que sera comparado con la variable "P" que nos indica la probabilidad
    de que exista un enlace entre dos nodos cualquiera. Se consideta un grafo aleatorio.
    """

    def ER(self, n, p, output_path=None):
        G = nx.Graph()

        for i in range(0, n):
            G.add_node(i)

        for i in range(0, n):
            for j in range(i + 1, n):
                if not G.has_node([i, j]) and random.uniform(0, 1) <= p:
                    G.add_edge(i, j)

        if output_path is not None:
            plt.title('Erdos-Renyi G = (n, p)  -->  G = ({}, {})'.format(n, p))
            nx.draw_spring(G, with_labels=True)
            plt.savefig(os.path.join(output_path, 'ER' + '_' + str(p) + '.png'))
            plt.clf()

            # exporting to geffy
            nx.write_gexf(G, os.path.join(output_path, 'ER' + '_' + str(p) + '.gexf'))

            adjacency_matrix = nx.to_numpy_matrix(G)
            max_k = n
            binomial_distribution = self.getBinomialDistribution(max_k, n, p)

            plotter = Plotter()
            plotter.histLinearScalePDF_ER(binomial_distribution, p, adjacency_matrix, output_path)

            # comienza desde 0 bin
            #plotter.histLinearScalePDF_ER_HISTO(binomial_distribution, p, adjacency_matrix, output_path)

        return G

    def getBinomialDistribution(self,k_max, n, p):
        binomial_dist = []

        for k in range(0, k_max):
            combinatorial_number = math.factorial(n - 1) / (math.factorial(int(k)) * math.factorial((n - 1) - int(k)))
            binomial_dist.append(combinatorial_number * (pow(p, k)) * pow(1 - p, n - 1 - k))

        return binomial_dist

if __name__ == '__main__':
    generator = ER()

    #ER: different values of "K" for G(N,K), or of "p" for G(N,p), e.g. p=0.00, 0.01, 0.02, 0.03, 0.05, 0.1
    generator.ER(50, 0.1, 'YOUR_OUTPUT_HERE')