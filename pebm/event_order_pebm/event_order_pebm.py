# Author: Christopher S. Parker (christopher.parker@ucl.ac.uk)
# P-EBM paper: https://www.biorxiv.org/content/10.1101/2022.07.10.499471v1.abstract
# Functionality inherited from KDE-EBM: https://github.com/ucl-pond/kde_ebm (Nick Firth, Neil Oxtoby, et al)

# License: TBC
import numpy as np
from random import sample, randrange
from itertools import chain
from numpy import ones, nan, prod, concatenate, log


class EventOrder_pebm(object):
    def __init__(self, ordering=None, n_biomarkers=None, score=None):
        super(EventOrder_pebm, self).__init__()
        if ordering is None and n_biomarkers is None:
            raise ValueError('EventOrder_pebm __init__ takes one argument,'
                             ' zero given')
        if ordering is None:
            self.ordering = [[i] for i in sample(range(n_biomarkers), n_biomarkers)]
            self.n_biomarkers = n_biomarkers
        else:
            self.ordering = ordering
            self.n_biomarkers = len(list(chain(*ordering)))
        self.score = score

    def score_ordering(self, prob_mat):
        event_order = self.ordering
        k = len(event_order) + 1
        p_perm = self.calc_perm_matrix(prob_mat)

        likelihood = np.sum(np.log(np.sum((1. / k) * p_perm, 1)))
        self.score = likelihood

        return likelihood

    def calc_perm_matrix(self, prob_mat):
        pxe = prob_mat[:, :, 1]
        pxne = prob_mat[:, :, 0]

        # pre-computing stages event/not-event indices
        event_order = self.ordering
        k = len(event_order) + 1
        pxes_inds_stages = [list(chain.from_iterable(event_order[:k_i])) for k_i in range(k)]
        pxnes_inds_stages = [list(chain.from_iterable(event_order[k_i:])) for k_i in range(k)]

        # compute permutation matrix
        p_perm = np.zeros((prob_mat.shape[0], k))
        for i in range(k):
            p_perm[:, i] = np.prod(pxe[:, pxes_inds_stages[i]], 1)*np.prod(pxne[:, pxnes_inds_stages[i]], 1) # each subjects (row) likelihood at stage k (column)
        return p_perm

    def stage_data(self, prob_mat):
        pxe = prob_mat[:, :, 1]
        pxne = prob_mat[:, :, 0]

        # pre-computing stages event/not-event indices
        event_order = self.ordering
        k = len(event_order) + 1
        pxes_inds_stages = [list(chain.from_iterable(event_order[:k_i])) for k_i in range(k)]
        pxnes_inds_stages = [list(chain.from_iterable(event_order[k_i:])) for k_i in range(k)]

        n_particp = prob_mat.shape[0]
        stage_likelihoods = np.empty((n_particp, k))
        for i in range(k):
            stage_likelihoods[:, i] = np.prod(pxe[:, pxes_inds_stages[i]], 1)*np.prod(pxne[:, pxnes_inds_stages[i]], 1)
        # Maximum-Likelihood stage
        stages = np.argmax(stage_likelihoods, axis=1)
        return stages, stage_likelihoods

    def swap_events(self):
        event_order = self.ordering
        # remove random node
        node_set = [item for sublist in event_order for item in sublist]
        n_nodes = len(node_set)
        random_node_id = randrange(0, n_nodes)
        new_event_order = [[ele for ele in sub if ele != random_node_id] for sub in event_order]
        new_event_order = [sub for sub in new_event_order if sub != []]
        # place at new position
        n_poss_pos = len(new_event_order) * 2 + 1
        random_pos = randrange(0, n_poss_pos)
        if random_pos % 2 == 0:
            random_insert_pos = int(random_pos / 2)
            new_event_order = new_event_order[:random_insert_pos] + [[random_node_id]] + new_event_order[
                                                                         random_insert_pos:]  # inserts *before* pos index
        elif random_pos % 2 == 1:
            random_add_pos = int((random_pos - 1) / 2)
            new_event_order[random_add_pos] += [random_node_id]
        return EventOrder_pebm(ordering=new_event_order)

    def __eq__(self, other):
        return self.ordering == other.ordering

    def __hash__(self):
        return hash(('ordering', str(self.ordering)))

    def __lt__(self, other):
        if self.score is None and other.score is None:
            raise ValueError('Cannot compare unscored orderings')
        if self.score < other.score:
            return True
        return False

    def __gt__(self, other):
        if self.score is None and other.score is None:
            raise ValueError('Cannot compare unscored orderings')
        if self.score > other.score:
            return True
        return False

    def __add__(self, other):
        if self.score is None and other.score is None:
            raise ValueError('Cannot add unscored orderings')
        return self.score + other.score

    def __sub__(self, other):
        if self.score is None and other.score is None:
            raise ValueError('Cannot subtract unscored orderings')
        return self.score - other.score

    def __repr__(self):
        return 'EventOrder(order=%r, score=%r)' % (self.ordering,
                                                   self.score)

    def __str__(self):
        return self.__repr__()

