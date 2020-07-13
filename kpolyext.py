# Thomas J. Savitsky
# July 10, 2020

import argparse
import sys
from igraph import *

class KPolyMatroid:
    """A k-polymatroid"""

    def __init__(self):
        self.maxn = 16  # matroids on 10 elements still haven't been enumerated
        self.flats = {}
        self.n = None
        self.flat_array = None
        self.flat_lookup = {} # the index of a mask in the flat array
        self.flat_containment_graph = None
        self.flat_cover_graph = None
        self.closure_memo = {}
        self.delta_memo = {}

    def __str__(self):
        ans = ''
        ary = self.get_flat_array()
        for f in ary:
            ans += format(f[1],'x') + ',' + str(f[0]) + ' '
        ans = ans.rstrip()
        return ans

    def get_flat_array(self):
        """an array of tuples of the form (rank,mask) where rank is the rank
            of the flat with bitmask mask.  ordered lexicographically.
            as a side effect, constructs self.flat_lookup"""
        if self.flat_array is not None:
            return self.flat_array
        ary = []
        for mask in self.flats:
            ary.append((self.flats[mask], mask))
        ary.sort()
        self.flat_array = ary

        for i in range(len(ary)):
            self.flat_lookup[ary[i][1]] = i

        return ary

    def get_flat_containment_graph(self):
        """a digraph with vertex-set equal to indices of the flat array.
            edge (i,j) exists iff the flat with index i is a proper subset
            of the flat with index j"""
        if self.flat_containment_graph is not None:
            return self.flat_containment_graph
        ary = self.get_flat_array()
        numflats = len(ary)
        g = Graph(n=numflats, directed=True)
    
        for i in range(numflats):
            for j in range(numflats):
                if i==j:
                    continue
                if (ary[i][1] & ary[j][1]) == ary[i][1]:
                    g.add_edge(i,j)

        self.flat_containment_graph = g
        return g
        
    def get_flat_cover_graph(self):
        """a digraph with vertex-set equal to indices of the flat array.
            edge (i,j) exists iff the flat with index i is covered by
            of the flat with index j in the lattice of flats"""
        if self.flat_cover_graph is not None:
            return self.flat_cover_graph
        ary = self.get_flat_array()
        cont = self.get_flat_containment_graph()
        numflats = cont.vcount()
        g = Graph(n=numflats, directed=True)
        for i in range(numflats):
            for j in cont.neighborhood(i, mode='out', mindist=1):
                iscovered = True
                for k in cont.neighborhood(i, mode='out', mindist=1):
                    if j == k:
                        continue
                    if (ary[k][1] & ary[j][1]) == ary[k][1]:
                        iscovered = False
                        break
                if iscovered:
                    g.add_edge(i,j)
        self.flat_cover_graph = g
        return g

    def read_flats_str(self, flats_str):
        """Format is a sequence of bitmasks of flats and ranks
           For example:
                0,0 1,2 2,2 3,3
            is two lines placed in a plane"""
        maxrank = -1
        fs = flats_str.split()
        for f in fs:
            flat = f.split(',')
            bitmask = int('0x'+flat[0], 16)
            rank = int(flat[1])
            kpm.flats[bitmask] = rank
            if rank > maxrank:
                rank = maxrank
                self.n = bin(bitmask).count("1")

    def relabel_flat(self, mask, label):
        ans = 0
        for i in range(self.n):
            if (1<<i) & mask:
                ans |= 1 << label[i]
        return ans

    def canonical_label(self):
        """Returns a copy of self that is canonically labeled."""
        ary = self.get_flat_array()

        # create a bipartite graph representing the polymatroid
        g = Graph()
        g.add_vertices(self.n + len(ary))
        colors = [0] * self.n
        for i in range(len(ary)):
            r, mask = ary[i]
            colors.append(r+1) #color flats by their rank + 1
            for j in range(self.n):
                if (1<<j) & mask:
                    g.add_edge(j,self.n+i)

        labeling = g.canonical_permutation(color=colors)

        canon = KPolyMatroid()
        canon.n = self.n
        for f in self.flats:
            canon.flats[self.relabel_flat(f, labeling)] = self.flats[f]
        return canon

    def closure(self, mask):
        if mask in self.flats:
            return mask
        if mask in self.closure_memo:
            return self.closure_memo[mask]
        closure = (1 << self.maxn) - 1
        for f in self.flats:
            if (mask & f == mask):
                closure &= f
        self.closure_memo[mask] = closure
        return closure

    def delta(self,f,g):
        """calculate modular defect"""
        if (f,g) in self.delta_memo:
            return self.delta_memo[(f,g)]
        md = (self.flats[f] + self.flats[g] - self.flats[f&g] -
                self.flats[self.closure(f|g)])
        self.delta_memo[(f,g)] = md
        self.delta_memo[(g,f)] = md
        return md

    def canonical_deletion(self):
        """returns the polymatroid obtained by:
            deleting element n-1 and canonically labelling this deletion"""
        deletion = KPolyMatroid()
        deletion.n = self.n-1

        for f in self.flats:
            if f & (1<<self.n-1):
                df = f & ~(1<<self.n-1)
                if df in self.flats:
                    continue
                else:
                    deletion.flats[df] = self.flats[f]
            else:
                deletion.flats[f] = self.flats[f]
            
        return deletion.canonical_label()

    def check_mu(self, mu):
        """returns True if the dictionary mu leads to a single-element
            extension"""
        ary = self.get_flat_array()
        numflats = len(ary)

        # condition I of Th. 3.3 could fail
        for i in range(numflats):
            frank, fmask = ary[i]
            for j in range(i+1, numflats):
                grank, gmask = ary[j]
                intmask = fmask & gmask
                intrank = self.flats[intmask]
                unionmask = self.closure(fmask|gmask)
                unionrank = self.flats[unionmask]
                if (mu[intmask] + mu[unionmask] - self.delta(fmask,gmask) >
                    mu[fmask] + mu[gmask]):
                    return False
        return True

    def produce_ext_from_mu(self, mu):
        ext = KPolyMatroid()
        ext.n = self.n + 1
        e = 1<<self.n
        ary = self.get_flat_array()
        covers = self.get_flat_cover_graph()

        for f in mu:
            frank = self.flats[f]
            if mu[f] == 0:
                ext.flats[f|e] = frank
                continue

            ext.flats[f] = frank

            i = self.flat_lookup[f] # an index into the flat array
            addflat = True
            for j in covers.neighborhood(i, mode='out', mindist=1):
                grank, gmask = ary[j]
                if frank + mu[f] == grank + mu[gmask]:
                    addflat = False
                    break
            if addflat:
                ext.flats[f|e] = frank + mu[f]
                
        return ext

    def force_mu(self, c, mu):
        """certain flats are required to be in mu_c"""
        ary = self.get_flat_array()
        cont = self.get_flat_containment_graph()
        covers = self.get_flat_cover_graph()

        # ensure condition II of Th. 3.3
        for i in range(len(ary)):
            frank, fmask = ary[i]
            if fmask in mu:
                continue
            for j in covers.neighborhood(i, mode='out', mindist=1):
                grank, gmask = ary[j]
                if gmask not in mu:
                    continue
                if (c == grank + mu[gmask] - frank):
                    mu[fmask] = c
                    break

        # ensure condition I of Th. 3.3
        forced = 1
        while (forced > 0):
            forced = 0
            mucpy = mu.copy()
            for f in mucpy:
                for g in mucpy:
                    if (f >= g) or (f&g in mu):
                        continue
                    if c == (mu[f] + mu[g] + self.delta(f,g) -
                              mu[self.closure(f|g)]):
                        forced += 1
                        mu[f&g] = c
                        # ensure condition III
                        i = self.flat_lookup[f&g]
                        for j in cont.neighborhood(i, mode='out', mindist=1):
                            grank, gmask = ary[j]
                            if gmask not in mu:
                                mu[gmask] = c

        return mu

    def get_mu_graph(self, c, mu):
        """The independent sets in this graph are candidates
            for the minimal elements of mu_c."""
        ary = self.get_flat_array()
        numflats = len(ary)
        g = Graph()

        for i in range(numflats):
            g.add_vertex(i)

        deletions = {}
        for i in range(numflats):
            frank, fmask = ary[i]
            if fmask in mu:
                deletions[i] = True

        for i in range(numflats):
            frank, fmask = ary[i]
            if fmask in mu:
                for j in range(numflats):
                    if (i==j) or (j in deletions):
                        continue
                    grank, gmask = ary[j]
                    k = (mu[fmask] + c + self.delta(fmask,gmask) -
                            mu[self.closure(fmask|gmask)])
                    if (k < c):
                        # g cannot be in mu_c
                        deletions[j] = True
                    elif (k == c) and (fmask & gmask != gmask):
                        # g cannot be a minimal member of mu_c
                        deletions[j] = True
            else:
                for j in range(i+1, numflats):
                    if (i in deletions) or (j in deletions):
                        continue
                    grank, gmask = ary[j]                    
                    if ((fmask & gmask == gmask) or (fmask & gmask == fmask)):
                        g.add_edge(i,j)
                        continue
                    union = self.closure(fmask|gmask)
                    if union in mu:
                        mu_union = mu[union]
                    else:
                        mu_union = c
                    if 2*c + self.delta(fmask,gmask) - mu_union <= c:
                        g.add_edge(i,j)

        g.delete_vertices(deletions.keys())
        return g    

    def ext_generator(self, clevel, maxc, mu):
        """recursively yield candidates for mu"""
        if clevel == maxc:
            for f in self.flats:
                if f not in mu:
                    mu[f] = clevel
            yield mu
        else:
            ary = self.get_flat_array()
            cont = self.get_flat_containment_graph()
            # certain flats will be forced to have mu[f] = clevel
            mu = self.force_mu(clevel, mu)
            g = self.get_mu_graph(clevel, mu)
            I = g.independent_vertex_sets()
            I.append(())  # add the empty set
            for X in I:
                minimals = g.vs(X)['name']    
                nextmu = mu.copy()
                for i in minimals:
                    frank, fmask = ary[i]
                    nextmu[fmask] = clevel
                    for j in cont.neighborhood(i, mode='out', mindist=1):
                        grank, gmask = ary[j]
                        if gmask not in nextmu:
                            nextmu[gmask] = clevel
                yield from self.ext_generator(clevel+1, maxc, nextmu)

    def extend(self, all, c, outf):
        """find all single-element extensions"""
        collection = {}
        for mu in self.ext_generator(0, c, {}):
            if self.check_mu(mu):
                ext = self.produce_ext_from_mu(mu)
                if all:
                    outf.write(str(ext)+'\n')
                else:
                    canonext = ext.canonical_label()
                    canondel = canonext.canonical_deletion()
                    if (str(self) == str(canondel) and
                        str(canonext) not in collection):
                            collection[str(canonext)] = True
                            outf.write(str(canonext)+'\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='kpolyext.py',
        description='Generate single-element extensions of a k-polymatroid')
    parser.add_argument('--label', help='Just canonically label the '\
        'polymatroids.  Do not generate extensions.', action='store_true',
        default=False)
    parser.add_argument('--all', help='Generate all extensions without '\
        'doing isomorphism checks. Defaults to False.', action='store_true',
        default=False)
    parser.add_argument('-c', help='The maximum rank of the new element. '\
            'Defaults to 1.', action='store',
            dest='max_rank', default=1, type=int)
    parser.add_argument('--version', action='version', version='%(prog)s 0.01')
    parser.add_argument('infile', help='name of input file',
        nargs='?', default=None, action='store')
    parser.add_argument('outfile', help='name of output file',
        nargs='?', default=None, action='store')

    results = parser.parse_args()

    if results.infile is None:
        inf = sys.stdin
    else:
        inf = open(results.infile, 'r')

    if results.outfile is None:
        outf = sys.stdout
    else:
        outf = open(results.outfile, 'w')

    for line in inf.readlines():
        kpm = KPolyMatroid()
        kpm.read_flats_str(line.rstrip())
        if results.label:
            canon = kpm.canonical_label()
            outf.write(str(canon)+'\n')
            continue
        kpm.extend(results.all, results.max_rank, outf)

    inf.close()
    outf.close()
