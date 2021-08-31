import numpy as np

class trellis(object):
    def __init__(self, m, g1, g0):
        self.numInputSymbols = pow(2, len(g1))
        self.numOutputSymbols = len(g1) + 1
        self.numStates = pow(2, m - 1)
        self.cl = m
        states = np.arange(self.numStates)
        l = np.arange(self.numInputSymbols)
        self.g0 = int(str(g0), 8)
        self.g1 = [int(str(elem), 8) for elem in g1]
        self.nextStates = np.ones((self.numStates, self.numInputSymbols))
        self.outputs = np.ones((self.numStates, self.numInputSymbols))
        for elemstate in states:
            for eleml in l:
                lastst = (sum([
                    int(i) for i in bin(self.g0 & (elemstate ^ (eleml <<
                                                                (m - 1))))[2:]
                ]) % 2 << (m - 1)) ^ elemstate
                newst = lastst >> 1
                self.nextStates[elemstate][eleml] = newst
                first = eleml << len(g1)
                for i in range(len(g1)):
                    first = first ^ ((
                        sum([int(j)
                             for j in bin(self.g1[i] & lastst)[2:]]) % 2) << i)
                self.outputs[elemstate][eleml] = first
                
def find_minweight_path(trellis, StartState, FinishState, T, maxZerosStates,
                        nZerosState):
    ## we want to find path with minimum weight of length T
    
    hdout = [0, 1, 1, 2]
    nStates = trellis.numStates
    wmat = 4 * T * np.ones((nStates, T + 2)) 
    wmat[FinishState, T + 1] = 0
    smat = -1 * np.ones((nStates, T + 2))
    smat[FinishState, T + 1] = FinishState
    chdist = 4 * T * np.ones((nStates, T + 1))
    curstates = -1 * np.ones(nStates)
    curstates[FinishState] = FinishState
    trmat = []
    hdmat = []

    for ii in np.arange(T, 0, -1):
        cwmat = 4 * T * np.ones((nStates, nStates))
        csmat = -1 * np.ones((nStates, nStates))

        pos = np.where(curstates != -1)
        curst = curstates[pos]

        for eljj in curst:
            curState = int(eljj)
            ind = trellis.nextStates == curState
            crosind = np.bitwise_and(ind[:, 0], ind[:, 1])
            ind0 = np.bitwise_xor(ind[:, 0], crosind)
            ind1 = np.bitwise_xor(ind[:, 1], crosind)
            if (sum(ind0.astype(int)) != 0):
                pos0 = np.where(ind0)
                csmat[pos0, curState] = [pos0[0]]
                w0 = wmat[curState][ii + 1]
                hdc0 = 4 * T * np.ones((nStates, 1))
                hdc0[pos0] = hdout[int(trellis.outputs[pos0, 0])]
                hdc0[pos0] = hdc0[pos0] + w0                    
                cwmat[:, curState] = hdc0[:, 0]
            if (sum(ind1.astype(int)) != 0):
                pos1 = np.where(ind1)
                csmat[pos1, curState] = [pos1[0]]
                w1 = wmat[curState][ii + 1]
                hdc1 = 4 * T * np.ones((nStates, 1))
                hdc1[pos1] = hdout[int(trellis.outputs[pos1, 1])]
                hdc01 = cwmat[:, curState]
                hdc01[pos1] = hdc1[pos1] + w1
                cwmat[:, curState] = hdc01
            if(curState==-1):
                print('Fucking chit!!!!!!!!!!!')
        if (csmat[0][0] == 0 and T != 1):
            if (nZerosState != maxZerosStates):
                nZerosState = nZerosState + 1
            else:
                csmat[0][0] = -1
                cwmat[0][0] = 4 * T
                nZerosState = 0
        if (len(trmat)):
            trmat = np.hstack((csmat, trmat))
        else:
            trmat = csmat
        if (len(hdmat)):
            hdmat = np.hstack((cwmat, hdmat))
        else:
            hdmat = cwmat
        pos = csmat != -1
        xcol = pos[:, 0]
        ocol = xcol
        for kk in range(1, nStates):
            xcol = np.bitwise_xor(xcol, pos[:, kk])
            ocol = np.bitwise_or(ocol, pos[:, kk])
        npos = np.sum(pos, axis=0)
        cs = -1 * np.ones(nStates)
        if (sum(np.bitwise_xor(xcol, ocol)) == 0):
            posc = np.where(npos > 0)[0]
            cw = 4 * T * np.ones(nStates)
            for elkk in posc:
                pos1 = np.where(pos[:, elkk])
                cs[pos1[0]] = csmat[pos1[0], elkk]
                cw[pos1[0]] = cwmat[pos1[0], elkk]
        else:
            wmin = np.min(cwmat, axis=1)
            cs[np.where(ocol)] = np.where(ocol)[0]
            cw = wmin
        if (ii == 1):
            if (cs[StartState] == StartState):
                cs = -1 * np.ones(nStates)
                cs[StartState] = StartState
                tmp = cw[StartState]
                cw = 4 * T * np.ones(nStates)
                cw[StartState] = tmp
            else:
                MinWeight = -1
                cs = -1 * np.ones(nStates)
                cw = 4 * T * np.ones(nStates)
                return [MinWeight, trmat, hdmat]
        smat[:, ii] = cs
        wmat[:, ii] = cw
        curstates = smat[:, ii]
    MinWeight = wmat[StartState, 1]
    if (MinWeight == 4 * T):
        MinWeight = -1
    return [MinWeight, trmat, hdmat, smat, wmat]

def find_spectrum(trellis, StartState, FinishState, T):
    ## we want to find spectrum with minimum weight of length T
    
    hdout = [0, 1, 1, 2]
    nStates = trellis.numStates
    nextst = trellis.nextStates.astype(int)
    outs = trellis.outputs.astype(int)
    nextst[0][0] = -1
    
    dwght = {}
    dcounts = {}
    
    states = []
    weights = []
    ind0 = nextst[:,0] == FinishState
    ind1 = nextst[:,1] == FinishState
    if(sum(ind0.astype(int)) != 0):
        ind = np.where(ind0)[0][0]
        states.append(ind)
        weights.append(hdout[outs[:,0][ind]])
    if(sum(ind1.astype(int)) != 0):
        ind = np.where(ind1)[0][0]
        states.append(ind)
        weights.append(hdout[outs[:,1][ind]])
    for i in range(T-1):
        cstates = states.copy()
        cweights = weights.copy()
        states = []
        weights = []
        for j, elem in enumerate(cstates):
            ind0 = nextst[:, 0] == elem
            w = cweights[j]
            if(sum(ind0.astype(int)) != 0):
                ind = np.where(ind0)[0][0]
                states.append(ind)
                weights.append(w+hdout[outs[:,0][ind]])
            ind1 = nextst[:, 1] == elem
            if(sum(ind1.astype(int)) != 0):
                ind = np.where(ind1)[0][0]
                states.append(ind)
                weights.append(w+hdout[outs[:,1][ind]])
        difw = []
        for j,elem in enumerate(states):
            if(elem == StartState):
                difw.append(weights[j])
        if(difw):
            difw = np.array(difw)
            uniq, counts = np.unique(difw, return_counts=True)
            dwght[i+2] = uniq
            dcounts[i+2] = counts
    return dwght, dcounts

def find_minweight_cword(trellis, FromState, T, trmat, hdmat, hd):
    hdout = np.array([0, 1, 1, 2])
    nStates = trellis.numStates
    tuplets = -1 * np.ones(T)
    infbits = -1 * np.ones(T)
    states = -1 * np.ones(T)
    hd1 = 0
    if (hd >= 0):
        for ii in range(1, T + 1):
            csmat = trmat[:, nStates * (ii - 1) + 0:nStates * ii]
            cwmat = hdmat[:, nStates * (ii - 1) + 0:nStates * ii]
            ind = csmat == FromState
            ind = ind[FromState]
            if (np.sum(ind) == 1):
                pos = np.where(ind)
                ToState = pos[0][0]
            else:
                pos = np.where(ind)
                wht = cwmat[FromState, :]
                if (wht[pos[0][0]] > wht[pos[0][1]]):
                    ToState = pos[0][1]
                else:
                    ToState = pos[0][0]
            states[ii - 1] = ToState
            ind = trellis.nextStates == ToState
            ind = ind[FromState, :]
            pos = np.where(ind)[0]
            if (pos[0] == 1):
                inbit = 1
            else:
                inbit = 0
            ctuplet = trellis.outputs[FromState, inbit]
            tuplets[ii - 1] = ctuplet
            infbits[ii - 1] = inbit
            hd1 = hd1 + hdout[int(ctuplet)]
            FromState = ToState
        if (hd != hd1):
            hd1 = -1
    else:
        hd1 = hd
    return [hd1, infbits, tuplets]