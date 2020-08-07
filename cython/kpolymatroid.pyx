#!python
#cython: language_level=3

# Thomas J. Savitsky
# August 7, 2020

from igraph import *

class KPolyMatroid:
    """A k-polymatroid"""

    def __init__(self):
        self.flats = {}
        self.n = None
        self.numflats = None
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
        self.numflats = len(ary)
        for i in range(self.numflats):
            self.flat_lookup[ary[i][1]] = i
        return ary

    def get_flat_containment_graph(self):
        """a digraph with vertex-set equal to indices of the flat array.
            edge (i,j) exists iff the flat with index i is a proper subset
            of the flat with index j"""
        if self.flat_containment_graph is not None:
            return self.flat_containment_graph
        ary = self.get_flat_array()
        numflats = self.numflats
        g = Graph(n=numflats, directed=True)
        edges = []
    
        for i in range(numflats):
            for j in range(numflats):
                if i==j:
                    continue
                if (ary[i][1] & ary[j][1]) == ary[i][1]:
                    edges.append((i,j))

        g.add_edges(edges)
        self.flat_containment_graph = g
        return g
        
    def get_flat_cover_graph(self):
        """a digraph with vertex-set equal to indices of the flat array.
            edge (i,j) exists iff the flat with index i is covered by
            of the flat with index j in the lattice of flats"""
        if self.flat_cover_graph is not None:
            return self.flat_cover_graph
        ary = self.get_flat_array()
        numflats = self.numflats
        cont = self.get_flat_containment_graph()
        g = Graph(n=numflats, directed=True)
        edges = []
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
                    edges.append((i,j))
        g.add_edges(edges)
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
            self.flats[bitmask] = rank
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
        numflats = self.numflats

        # create a bipartite graph representing the polymatroid
        g = Graph()
        g.add_vertices(self.n + numflats)
        edges = []
        colors = [0] * self.n
        for i in range(numflats):
            r, mask = ary[i]
            colors.append(r+1) #color flats by their rank + 1
            for j in range(self.n):
                if (1<<j) & mask:
                    edges.append((j,self.n+i))

        g.add_edges(edges) #massive speedup compared to adding one-by-one
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
        closure = (1 << self.n) - 1
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

    def check_mu(self, mu, c):
        """returns True if the mu constructed so far
            could lead to a single-element extension"""
        # only condition I of Th. 3.3 could fail at this point, and
        # only for flats just added to mu_c
        mu_c = {}
        for f in mu:
            if mu[f] == c:
                mu_c[f] = c
    
        for f in mu:
            for g in mu_c:
                if f == g:
                    continue
                k = mu[f] + c - mu[self.closure(f|g)] + self.delta(f,g)
                if (f&g) in mu_c:
                    if c > k:
                        return False
                elif c+1 > k:  #if c <= k, mu[f&g] can be assigned later
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
        numflats = self.numflats
        cont = self.get_flat_containment_graph()
        covers = self.get_flat_cover_graph()

        mucpy = mu.copy()
        # ensure condition II of Th. 3.3
        for i in range(numflats):
            frank, fmask = ary[i]
            if fmask in mucpy:
                continue
            for j in covers.neighborhood(i, mode='out', mindist=1):
                grank, gmask = ary[j]
                if gmask not in mucpy:
                    continue
                if (c == grank + mucpy[gmask] - frank):
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
                    if c == (mucpy[f] + mucpy[g] + self.delta(f,g) -
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
        """The independent sets in this graph are candidates for
            minimal elements of mu_c that have not already been forced in."""
        ary = self.get_flat_array()
        numflats = self.numflats
        edges = []
        gph = Graph()

        gph.add_vertices(range(numflats)) # we need the vertex labels

        deletions = {}
        for i in range(numflats):
            frank, fmask = ary[i]
            if fmask in mu:
                deletions[i] = True

        for i in range(numflats):
            frank, fmask = ary[i]
            if fmask in mu:
                for j in range(numflats):
                    if j in deletions:
                        continue
                    grank, gmask = ary[j]
                    k = (mu[fmask] + c + self.delta(fmask,gmask) -
                            mu[self.closure(fmask|gmask)])
                    if (c > k):
                        # f,g fails condition I if g is in mu_c
                        deletions[j] = True
                    elif (c == k) and (fmask & gmask != gmask):
                        # f,g fails if g is a minimal member of mu_c
                        deletions[j] = True

        for i in range(numflats):
            if i in deletions:
                continue
            frank, fmask = ary[i]
            for j in range(i+1, numflats):
                if j in deletions:
                    continue
                grank, gmask = ary[j]                    
                if ((fmask & gmask == gmask) or (fmask & gmask == fmask)):
                    edges.append((i,j))
                    continue
                union = self.closure(fmask|gmask)
                if union in mu:
                    mu_union = mu[union]
                else:
                    mu_union = c
                if c+1 > 2*c + self.delta(fmask,gmask) - mu_union:
                    edges.append((i,j))

        gph.add_edges(edges)
        gph.delete_vertices(deletions.keys())
        return gph

    def ext_generator(self, c, maxc, mu):
        """recursively yield candidates for mu"""
        if c == maxc:
            for f in self.flats:
                if f not in mu:
                    mu[f] = c
            yield mu
        else:
            ary = self.get_flat_array()
            cont = self.get_flat_containment_graph()
            # certain flats will be forced to have mu[f] = c
            mu = self.force_mu(c, mu)
            g = self.get_mu_graph(c, mu)
            I = g.independent_vertex_sets()
            I.append(())  # add the empty set
            for X in I:
                minimals = g.vs(X)['name']    
                nextmu = mu.copy()
                for i in minimals:
                    frank, fmask = ary[i]
                    nextmu[fmask] = c
                    for j in cont.neighborhood(i, mode='out', mindist=1):
                        grank, gmask = ary[j]
                        if gmask not in nextmu:
                            nextmu[gmask] = c
                if self.check_mu(nextmu, c):
                    yield from self.ext_generator(c+1, maxc, nextmu)

    def check_mu_debug(self, mu):
        """returns True if the dictionary mu leads to a single-element
            extension"""
        ary = self.get_flat_array()
        numflats = self.numflats
        for i in range(numflats):
            frank, fmask = ary[i]
            for j in range(numflats):
                grank, gmask = ary[j]
                # condition I
                intmask = fmask & gmask
                intrank = self.flats[intmask]
                unionmask = self.closure(fmask|gmask)
                unionrank = self.flats[unionmask]
                if (mu[intmask] + mu[unionmask] - self.delta(fmask,gmask) >
                    mu[fmask] + mu[gmask]):
                    return False
                # condition II
                if (fmask & gmask == fmask):
                    if (self.flats[fmask] + mu[fmask] > self.flats[gmask]
                        + mu[gmask]):
                        return False
                # condition III
                if (fmask & gmask == fmask):
                    if mu[gmask] > mu[fmask]:
                        return False
        return True

    def extend(self, all, c, max_r, outf):
        """find all single-element extensions"""
        collection = {}
        for mu in self.ext_generator(0, c, {}):
            #if not self.check_mu_debug(mu):
            #    print("Bad extension.  c=", c, "mu=", mu)
            #    exit(1)
            ext = self.produce_ext_from_mu(mu)
            if (max_r is not None):
                if ext.flats[(1<<ext.n)-1] > max_r:
                    continue
            if all:
                outf.write(str(ext)+'\n')
            else:
                canonext = ext.canonical_label()
                canondel = canonext.canonical_deletion()
                if (str(self) == str(canondel) and
                        str(canonext) not in collection):
                    collection[str(canonext)] = True
                    outf.write(str(canonext)+'\n')
