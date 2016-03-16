import sys
import random
import time
import oph

random.seed()
debug = True if 'd' in sys.argv else False
f = open( 'test instances/tsiligirides_problem_3_budget_015.txt' )
tmax = int( f.readline().split()[0] )
tmax = 25
cpoints = []
for i in xrange( 33 ):
    cpoints.append( tuple( [ float( x ) for x in f.readline().split() ] + [ i, 0 ] ) )
start_point = cpoints.pop( 0 )
end_point = cpoints.pop( 0 )
popsize = 200
genlimit = 100
kt = 5
isigma = 15
musigma = 7
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

#fitness will take a set of weights and return a tuple containing the fitness and the best path
def fitness( chrom ):
    augs = []
    for i in xrange( len( cpoints ) ):
        augs.append( ( cpoints[ i ][0],
                       cpoints[ i ][1],
                       cpoints[ i ][2],
                       cpoints[ i ][3], 
                       cpoints[ i ][4] + chrom[ i ][4] ) )
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
        print len( cpoints )
        print [ cpoints[ x[3] - 2 ] for x in best[ 1:len( best ) - 1 ] ]
        print [ cpoints[ x[3] - 2 ][2] for x in best[ 1:len( best ) - 1 ] ]
        print ( sum( [ cpoints[ x[3] - 2 ][2] for x in best[ 1:len( best ) - 1 ] ] ), best )
    return ( sum( [ cpoints[ x[3] - 2 ][2] for x in best[ 1:len( best ) - 1 ] ] ), best )

def crossover( c1, c2 ):
    assert( len( c1 ) == len( c2 ) )
    point = random.randrange( len( c1 ) )
    first = random.randrange( 2 )
    if( first ):
        return c1[:point] + c2[point:]
    else:
        return c2[:point] + c1[point:]

def mutate( chrom ):
    for i in xrange( len( chrom ) ):
        mulist = []
        for j in xrange( 4, len( chrom[ i ] ) ):
            if( random.randrange( mchance ) == 0 ):
                mulist.append( chrom[ i ][ j ] + random.gauss( 0, musigma ))
            else:
                mulist.append( chrom[ i ][ j ] )
        chrom[ i ] = tuple( [0, 0, 0, 0] + mulist )
    return chrom

start_time = time.clock()
#generate initial random population
pop = []
for i in xrange( popsize + elitismn ):
    chrom = []
    for j in xrange( len( cpoints ) ):
        chrom.append( ( 0, 0, 0, 0, random.gauss( 0, isigma ) ) )
    chrom = ( fitness( chrom )[0], chrom )
    while( i - j > 0 and j < elitismn and chrom > pop[ i - 1 - j ] ):
        j += 1
    pop.insert( i - j, chrom )
    #pop.append( ( fitness( chrom )[0], chrom ) )

bestfit = 0
for i in xrange( genlimit ):
    nextgen = []
    for j in xrange( popsize ):
        #select parents in k tournaments
        parents = sorted( random.sample( pop, kt ) )[ kt - 2: ] #optimize later
        #crossover and mutate
        offspring = mutate( crossover( parents[0][1], parents[1][1] ) )
        offspring = ( fitness( offspring )[0], offspring )
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
best_path = fitness( bestchrom[1] )[1]
print [ x[3] for x in best_path ]

print 'their shit:'
shit = oph.initialize( oph.ell_sub( tmax, start_point, end_point, cpoints )
, start_point, end_point, tmax )[0]
print 'fitness:', sum( [ x[2] for x in shit ] )
print 'my shit:'
shit2 = oph.ellinit_replacement( cpoints, start_point, end_point, tmax )
print 'fitness:', sum( [ x[2] for x in shit2 ] )
print 'checking correctness...',
total_distance = ( oph.distance( start_point, cpoints[ best_path[ 1                    ][3] - 2 ] ) + 
                   oph.distance( end_point,   cpoints[ best_path[ len( best_path ) - 2 ][3] - 2 ] ) )
for i in xrange( 1, len( best_path ) - 3 ):
    total_distance += oph.distance( cpoints[ best_path[ i     ][3] - 2 ], 
                                    cpoints[ best_path[ i + 1 ][3] - 2 ] )
print 'OK' if total_distance <= tmax else 'not OK'
print 'tmax:          ', tmax
print 'total distance:', total_distance
