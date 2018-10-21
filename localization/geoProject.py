from geoInterface import *
from geometry import *
from methods import *


class Project:
    def __init__(self, mode='2D', solver='LSE', detail=False):
        self.mode = mode
        self.solver = solver
        self.detail = detail
        self.AnchorDic = {}
        self.TargetDic = {}
        self.nt = 0

    def set_mode(self, mode):
        self.mode = mode

    def set_solver(self, sol):
        self.solver = sol

    def add_anchor(self, ID, loc):
        try:
            self.AnchorDic[ID]
            print(str(ID) + ':Anchor with same ID already exists')
            return
        except KeyError:
            a = Anchor(ID, point(loc))
            self.AnchorDic[ID] = a
        return a

    def add_target(self, ID=None):
        try:
            self.TargetDic[ID]
            print('Target with same ID already exists')
            return
        except:
            self.nt = self.nt + 1
            if ID:
                pass
            else:
                ID = 't' + str(self.nt)
            t = Target(ID)
            self.TargetDic[ID] = t
        return (t, ID)

    def solve(self, **kwargs):
        for tID in self.TargetDic.keys():
            tar = self.TargetDic[tID]
            cA = []
            for tup in tar.measures:
                landmark = tup[0]
                c = self.AnchorDic[landmark].loc
                d = tup[1]
                cA.append(circle(c, d))
            if self.solver == 'LSE':
                tar.loc = lse(cA, mode=self.mode, cons=False)
            elif self.solver == 'LSE_GC':
                try:
                    tar.loc = lse(cA, mode=self.mode, cons=True)
                except cornerCases as cc:
                    if cc.tag == 'Disjoint':
                        print(tar.ID + ' could not be localized by LSE_GC')
                    else:
                        print('Unknown Error in localizing ' + tar.ID)
            elif self.solver == 'CCA':
                if not self.detail:
                    tar.loc, n = CCA(cA, mode=self.mode, detail=False)
                    return n
                else:
                    tar.loc, n, P, iP = CCA(cA, mode=self.mode, detail=True)
                    return (n, P, iP)
