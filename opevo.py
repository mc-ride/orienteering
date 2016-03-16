import sys
import random
import time
import oph

#fitness will take a set s and a set of weights and return a tuple containing the fitness and the best path
def fitness( chrom, s, start_point, end_point ):
    augs = []
    for i in xrange( len( s ) ):
        augs.append( ( s[ i ][0],
                       s[ i ][1],
                       s[ i ][2],
                       s[ i ][3], 
                       s[ i ][4] + chrom[ i ] ) )
    if debug:
        print 'fitness---------------------------------'
        print 'augs:'
        print augs
    #best = oph.ellinit_replacement( augs, start_point, end_point, tmax )
    ellset = oph.ell_sub( tmax, start_point, end_point, augs )
    #best = oph.initialize( ellset, start_point, end_point, tmax )[0]
    best = oph.init_replacement( ellset, start_point, end_point, tmax )[0]
    if debug:
        print 'best:'
        print best
        print 'best real reward:'
        print [ x[3] for x in best ]
        print len( s )
        print [ s[ x[3] - 2 ] for x in best[ 1:len( best ) - 1 ] ]
        print [ s[ x[3] - 2 ][2] for x in best[ 1:len( best ) - 1 ] ]
        print ( sum( [ s[ x[3] - 2 ][2] for x in best[ 1:len( best ) - 1 ] ] ), best )
    return ( sum( [ s[ x[3] - 2 ][2] for x in best[ 1:len( best ) - 1 ] ] ), best )

def crossover( c1, c2 ):
    assert( len( c1 ) == len( c2 ) )
    point = random.randrange( len( c1 ) )
    first = random.randrange( 2 )
    if( first ):
        return c1[:point] + c2[point:]
    else:
        return c2[:point] + c1[point:]

def mutate( chrom, mchance, msigma ):
    return [ x + random.gauss( 0, msigma ) if random.randrange( mchance ) == 0  else 
             x for x in chrom ]

def run_alg( f, tmax, N ):
    random.seed()
    cpoints = []
    for i in xrange( N ):
        cpoints.append( tuple( [ float( x ) for x in f.readline().split() ] + [ i, 0 ] ) )
    start_point = cpoints.pop( 0 )
    end_point = cpoints.pop( 0 )
    popsize = 20
    genlimit = 50
    kt = 5
    isigma = 10
    msigma = 7
    mchance = 1
    elitismn = 2
    if( debug ):
        print 'data set size:', len( cpoints ) + 2
        print 'parameters:'
        print 'generations:     ', genlimit
        print 'population size: ', popsize
        print 'ktournament size:', kt
        print 'mutation chance: ', mchance
        print str( elitismn ) + '-elitism'

    start_time = time.clock()
    #generate initial random population
    pop = []
    for i in xrange( popsize + elitismn ):
        chrom = []
        for j in xrange( len( cpoints ) ):
            chrom.append( random.gauss( 0, isigma ) )
        chrom = ( fitness( chrom, cpoints, start_point, end_point )[0], chrom )
        while( i - j > 0 and j < elitismn and chrom > pop[ i - 1 - j ] ):
            j += 1
        pop.insert( i - j, chrom )

    bestfit = 0
    for i in xrange( genlimit ):
        nextgen = []
        for j in xrange( popsize ):
            #select parents in k tournaments
            parents = sorted( random.sample( pop, kt ) )[ kt - 2: ] #optimize later
            #crossover and mutate
            offspring = mutate( crossover( parents[0][1], parents[1][1] ), mchance, msigma )
            offspring = ( fitness( offspring, cpoints, start_point, end_point )[0], offspring )
            if( offspring[0] > bestfit ):
                bestfit = offspring[0]
                print bestfit
            if( elitismn > 0 and offspring > pop[ popsize ] ):
                l = 0
                while( l < elitismn and offspring > pop[ popsize + l ] ):
                    l += 1
                pop.insert( popsize + l, offspring )
                nextgen.append( pop.pop( popsize ) )
            else:
                nextgen.append( offspring )
        pop = nextgen + pop[ popsize: ]

    bestchrom = sorted( pop )[ popsize + elitismn - 1 ] #optimize later
    end_time = time.clock()

    print 'time:'
    print end_time - start_time
    print 'best fitness:'
    print bestchrom[0]
    print 'best path:'
    best_path = fitness( bestchrom[1], cpoints, start_point, end_point )[1]
    print [ x[3] for x in best_path ]

    print 'their stuff:'
    stuff = oph.initialize( oph.ell_sub( tmax, start_point, end_point, cpoints )
    , start_point, end_point, tmax )[0]
    print 'fitness:', sum( [ x[2] for x in stuff ] )
    print 'my stuff:'
    stuff2 = oph.ellinit_replacement( cpoints, start_point, end_point, tmax )
    print 'fitness:', sum( [ x[2] for x in stuff2 ] )
    print 'checking correctness...',
    total_distance = ( oph.distance( start_point, cpoints[ best_path[ 1                    ][3] - 2 ] ) + 
                       oph.distance( end_point,   cpoints[ best_path[ len( best_path ) - 2 ][3] - 2 ] ) )
    for i in xrange( 1, len( best_path ) - 3 ):
        total_distance += oph.distance( cpoints[ best_path[ i     ][3] - 2 ], 
                                        cpoints[ best_path[ i + 1 ][3] - 2 ] )
    print 'OK' if total_distance <= tmax else 'not OK'
    print 'tmax:          ', tmax
    print 'total distance:', total_distance

if( __name__ ==  '__main__' ):
    debug = True if 'd' in sys.argv else False
    #f = open( 'test instances/tsiligirides_problem_3_budget_015.txt' )
    f = open( sys.argv[1] )
    of = open( sys.argv[1] + '_results.dat', 'w' )
    #tmax = int( f.readline().split()[0] )
    tmax = int( sys.argv[2] )
    N = int( sys.argv[3] )
    run_alg( f, tmax, N )
