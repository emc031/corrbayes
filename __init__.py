
import sys
import re
import random
from math import isnan

import numpy as np
from numpy import array

from gvar import log,exp,evalcov,sqrt
from gvar.dataset import Dataset,avg_data
import gvar as gv

default_loosener=0.3
default_zerobuffer=0.1

def invertosc(c):
    """
    Swaps effect of oscillating and non-oscillating states in correlator "c" using the operation
    C(t) -> (-1)^t C(t).
    """
    return [ (-1)**(t+1) *c[t] for t in range(len(c)) ]
##


def superav(c,n):
    """
    Performs superaveraging operation on correlator "c" n times.
    The superaverage operator does this
     C(t) -> ( C(t)+C(t+1) ) / 2.
    """
    for m in range(n):
        c = [ (c[t]+c[t+1])/2 for t in range(len(c)-1)]
    return c
##


def superav2(c):
    """
    Asymmetric second-order superaverage of correlator "c";
     C(t) -> ( C(t-1) - 2C(t) + C(t+1) ) / 4 .
    """
    return [ (2*c[t]+c[t+1]+c[t+2])/4 for t in range(len(c)-2) ]
##


def effective_mass(c):
    """
    Returns effective mass of correlator "c"
     m_{eff}(t) = log( C(t) / C(t+1) ).
    """
    return [ log( sqrt( c[t]/c[t+1] )**2 ) for t in range(len(c)-1)]
##


def effective_amp(c): 
    """
    Returns effective amplitude
    a_{eff}(t) = sqrt( C(t) * exp( m_{eff}(t) * t ) )
    """
    m=effective_mass(c)
    return [ sqrt( c[t] * exp(m[t]*t) ) for t in range(len(c)-1) ]
##

def amp_superav(c): 
    """
    Takes superaverage of correlator c, then effective amplitude,
    then applies a correction factor sqrt( 2 / 1 + e^{-m_{eff}(t)} ).

    This is because taking the superaverage shifts the effective amplitude of a correlator,
    this correction cancels that shift.

    """
    c = superav(c)
    m=effective_mass(c)
    return [ sqrt( c[t] * exp(m[t]*t) * 2/(1+exp(-m[t])) )  for t in range(len(c)-1) ]
##

def amp_superav2(c): 
    """
    Takes superaverage of correlator c, then effective amplitude,
    then applies a correction factor sqrt( 4 / 1 + cosh(m_{eff}(t)) ).

    This is because taking the superaverage shifts the effective amplitude of a correlator,
    this correction cancels that shift.
    """
    c = superav2(c)
    m = effective_mass(c)
    return [ sqrt( c[t] * exp(m[t]*t) * 2/(1+(exp(-m[t])+exp(-2*m[t]))/2 ) ) for t in range(len(c)-2) ]
##


def ratio(C3, C2_1, C2_2):
    """
    Takes a 3-point correlator C3 and two 2-point correlators C2_1,C2_2, returns the ratio

    R = C_{3pt}(t)/ ( C_{2pt,1}(t) * C_{2pt,2}(T-t) )

    where T is the source/sink temporal separation.
    """
    T = C3.shape[0]
    return [C3[t]/(C2_1[T-t]*C2_2[t]) for t in range(T)]
##


def safelog(x, verbose=False):
    """
    Takes in gvar x, returns log(x). If x not suitable for log, returns log(1.0(9)).
    """
    logx = log(x)
    if isnan(logx.mean):
        if verbose: print('CorrBayes.safelog WARNING: invalid argument for log - replacing with log(1.0(9))')
        return log(gv.gvar(1.0,0.9))
    else: return logx
##


def dirtyfit(correlator,
             nexp,
             key,

             tcut=None,
             loosener=default_loosener,
             zero_buffer=default_zerobuffer,

             verbose=False
          ):
    """
    Takes in single correlator and associated key,
    returns dictionary of best guesses for correlator parameters.

    Requires data keys of the format "meson.ss", where "meson" labels the meson,
    and the "s"'s labels the source/sink combination,
    e.g. etac.ll is a local-local eta_c correlator.
    """

    if verbose: print('Performing dirty fit on correlator',key)

    Tlat = len(correlator)
    if not tcut: 
        tcut = int(Tlat/10.)
        if verbose: print('tcut set to ',tcut)

    if not loosener: loosener=default_loosener
    if not zero_buffer: zero_buffer=defalt_zerobuffer

    result = gv.BufferDict()

    # finding ground state energy
    mass = avg_data(
        effective_mass( superav2(correlator) )[tcut:int(Tlat/2)-tcut]
    ) * gv.gvar(1,loosener)
    if verbose: print('mass = ',mass)

    # finding ground state amplitude
    amp = avg_data(
        amp_superav2(correlator)[tcut:int(Tlat/2)-tcut]
    ) * gv.gvar(1,loosener)
    if verbose: print('amplitude = ',amp)

    # 'excited_correlator' = correlators with ground state term removed
    excited_correlator = [ correlator[t] - amp*exp(-mass*t) for t in range(Tlat) ]

    # finding first excited energy
    spectrum = np.mean(
        effective_mass( superav2(excited_correlator) )[ tcut : int(Tlat/2)-tcut ]
    ) * gv.gvar(1,loosener)

    spectrum = gv.gvar( zero_buffer/(1-spectrum.sdev/spectrum.mean), 
                        zero_buffer/(spectrum.mean/spectrum.sdev - 1) )
    # this transforms spectrum to a gvar that is 1 sigma
    # away from zero_buffer with the same fractional sdev

    # finding first excited state amplitude
    spectrum_amps = np.mean(
        amp_superav2( excited_correlator )[ tcut : int(Tlat/2)-tcut ]
    ) * gv.gvar(1,loosener)
    spectrum_amps = gv.gvar( zero_buffer/(1-spectrum_amps.sdev/spectrum_amps.mean), 
                        zero_buffer/(spectrum_amps.mean/spectrum_amps.sdev - 1) )
    # this transforms spectrum_amps to a gvar that is 1 sigma
    # away from zero_buffer with the same fractional sdev


    # building key
    try:
        meson = re.findall('^(.*)\.[a-zA-Z][a-zA-Z]$',key)[0] # M.ll -> M
        source = key.split('.')[-1][0] # M.ll -> l
    except IndexError:
        meson = key
        source=''
    if verbose: print('found meson label = ',meson,', source label = ',source)

    # building dictionary of results
    result.add('log'+meson+':a'+source,
               [ safelog(spectrum_amps) for i in range(nexp) ])
    result.add('logdE:'+meson,
               [ safelog(spectrum) for i in range(nexp) ])

    result['log'+meson+':a'+source][0] = safelog(amp)
    result['logdE:'+meson][0] = safelog(mass)

    # making guesses for oscillating states -
    # ground state amplitudes smaller and ground state energies larger than non-oscillating
    result.add('logo'+meson+':a'+source,
               [safelog( gv.gvar( amp.mean/2, amp.mean ) ) for i in range(nexp)] )
    result.add('logdE:o'+meson,
               [safelog(spectrum) for i in range(nexp)] )

    result['logdE:o'+meson][0] = safelog( gv.gvar( mass.mean*1.5, mass.mean ) )

    if verbose: print('result = ',result)

    return gv.add_parameter_parentheses(result)
##


def dirtyfit_3pt(cdict,
                 nexp,
                 key,
                 current,

                 tcut_3pt=None,
                 tcut_2pt=None,
                 loosener=default_loosener,

                 verbose=False
              ):
    """
    Takes in a dictionary of correlators, the key for a 3-point correlator and the name of the current,
    produces a dirty estimate of the 3-point transition amplitude by taking the ratio

    R = C_{3pt}(t)/ C_{2pt,1}(t) * C_{2pt,2}(T-t),

    where the C_{2pt,1/2} are the correlators for the two states on either side of the current.
    t is the timeslice of the current and T is the source/sink temporal separation.

    This ratio is equal to J/a1*a2, where J is the transition amplitude and a1/2 are the amplitudes
    of the 2-point correlators. Then we find J = R * a1 * a2.

    Requires 2-point data keys of the format "meson.ss", where "meson" labels the meson,
    and the "s"'s labels the source/sink combination,
    e.g. etac.ll is a local-local eta_c correlator.

    Requires 3-point data keys of the format "meson1.J.meson2_T{T}.ss".
    """

    if verbose: print('Performing dirty fit on correlator',key)

    if not loosener: loosener=default_loosener

    # fetching labels
    try:
        tag = re.findall('^(.*)_T\d\d?',key)[0] # M1.J.M2_T{T}.ll -> M1.J.M2
        meson1 = re.findall('^(.*)\.'+current,key)[0] # M1.J.M2_T{T}.ll -> M1
        meson2 = re.findall(current+'\.(.*)_T\d\d?',key)[0] # M1.J.M2_T{T}.ll -> M2

        source_ = re.findall('_T\d\d?\.[a-zA-Z]([a-zA-Z])$',key) # M1.J.M2_T{T}.ll -> l
        if len(source_)>0: source = source_[0]
        else: source=''

    except IndexError:
        print('The key',key,'is in the wrong format, it must look like M1.J.M2_T{T}.ll')
        sys.exit(1)

    T = int( re.findall('_T(\d\d?)',key)[0] ) # M1.J.M2_T{T}.ll -> {T}

    if verbose:
        print('found meson labels =',meson1,',',meson2,', sources =',source)
        print('found T =',T)

    c3 = cdict[key]

    if not tcut_3pt: 
        tcut_3pt = int(T/3.)
        if verbose: print('tcut_3pt set to',tcut_3pt)

    c2 = []; amp = []
    for meson in [meson1,meson2]:
        if verbose: print('finding amplitude for',meson,'...')

        # finding corresponding 2point correlators
        try:
            if source=='': c2.append( cdict[meson] )
            else: c2.append( cdict[meson+'.'+source+source] )
        except KeyError:
            print('cannot find correlator for',meson,'to go with',key)
            sys.exit(1)

        Tlat = len(c2[-1])

        if not tcut_2pt:
            tcut_2pt = int(Tlat/10.)
            if verbose: print('tcut_2pt set to',tcut_2pt)

        # finding ground state amplitude of 2pt correlators
        amp.append(
            avg_data(
                amp_superav2(c2[-1])[ tcut_2pt : int(Tlat/2)-tcut_2pt ]
            ) * gv.gvar(1,loosener)
        )

        if verbose: print('amp =',amp[-1])

    # J = approximation of 3-point transition amplitude
    J = avg_data(
        superav2( ratio(c3,c2[0],c2[1]) )[ tcut_3pt : T-tcut_3pt ]
    ) * np.product(amp) * gv.gvar(1,loosener)

    if verbose: print('J =',J)

    result = gv.BufferDict()
    for V in [ 'Vnn_', 'Von_', 'Vno_', 'Voo_' ]:
        result.add( V+tag,
                    [ [gv.gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)] )

    result['Vnn_'+tag][0][0] = J

    if verbose: print('result =',result)

    return gv.add_parameter_parentheses(result)
##



def get_prior(dset,
              Nsubset,
              nexp,
              currents=[],

              tcut_2pt=None,
              tcut_3pt=None,
              loosener=0.3,
              zero_buffer=0.1,

              verbose=False
          ):
    """
    Takes in a gvar.Dataset "dset", shaves off "Nsubset" data points,
    and uses those correlators to deduce priors for a fit of the correlators
    in dset.
    """

    assert type(dset)==gv.dataset.Dataset
    assert type(Nsubset)==int and Nsubset > 0
    assert type(nexp)==int and nexp > 0

    # pick out subset

    subset_index = random.sample( list(range(0,len(dset[list(dset.keys())[0]]))),
                                  Nsubset )
    if verbose: print('data point(s) selected for deducing prior:',subset_index)
    subset = { key :
               np.mean( [ dset[key][s] for s in subset_index ], axis=0 )
               for key in list(dset.keys()) } # averaged subset data

    # removing subset from dataset
    if verbose: print('removing these points from dataset...')
    for s in subset_index:
        for key in list(dset.keys()):
            dset[key] = np.concatenate( ( dset[key][:s], dset[key][s+1:] ) )


    # deduce priors from subset

    prior = gv.BufferDict()
    for key in list(subset.keys()):
        c = subset[key]

        # see if this key corresponds to a 3pt correlator, and save current name
        for current_i in currents:
            if re.search(current_i,key):
                this_current = current_i

        # do single correlator fit
        if 'this_current' in locals():
            cpriors = dirtyfit_3pt( cdict = subset,
                                    nexp = nexp,
                                    key = key,
                                    current = this_current,
                                    tcut_3pt = tcut_3pt,
                                    tcut_2pt = tcut_2pt,
                                    loosener = loosener,
                                    verbose = verbose )

            del this_current
        else:
            cpriors = dirtyfit( correlator = c,
                                nexp = nexp,
                                key = key,
                                tcut = tcut_2pt,
                                loosener = loosener,
                                zero_buffer = zero_buffer,
                                verbose = verbose )

        for key in list(cpriors.keys()): prior[key] = cpriors[key]

    return ( prior, dset )
##
