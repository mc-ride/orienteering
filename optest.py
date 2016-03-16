import time
import opevo

files = [ 'test instances/tsiligirides_problem_1_budget_15.txt',
          'test instances/tsiligirides_problem_2_budget_15.txt',
          'test instances/tsiligirides_problem_3_budget_015.txt',
          'test instances/set_66_1_015.txt',
          'test instances/set_64_1_15.txt' ]
tmaxs = [ range(  5,  85 + 1, 5 ),
          range( 15,  45 + 1, 5 ),
          range( 15, 110 + 1, 5 ),
          range(  5, 130 + 1, 5 ),
          range( 15,  80 + 1, 5 ) ]
Ns = [ 31, 21, 33, 66, 64 ]

test_runs = 30

assert( len( files ) == len( tmaxs ) and len( tmaxs ) == len( Ns ) )

for i in xrange( len( files ) ):
    f = open( files[ i ] )
    of = open( files[ i ][ :len( files[ i ] ) - 4 ] + '_results.dat', 'a' )

    of.write( time.asctime() + '\n' )
    of.write( 't avgfit avgtime bestfit\n' )
    for t in tmaxs[ i ]:
        fit_sum = float( 0 )
        time_sum = float( 0 )
        best_fit = 0
        for j in xrange( test_runs ):
            print 'TEST %i/%i' % ( j + 1, test_runs )
            f.seek( 0 ) 
            result = opevo.run_alg( f, t, Ns[ i ] )
            fit_sum += result[0]
            time_sum += result[1]
            best_fit = result[0] if result[0] > best_fit else best_fit
        #find avg fit, time, best fit then write to file
        of.write( ' '.join( [ str( x ) for x in [ t, fit_sum / test_runs, time_sum / test_runs,
            best_fit ] ] ) + '\n' )
    f.close()
    of.close()
