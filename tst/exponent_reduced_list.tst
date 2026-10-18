gap> START_TEST( "exponent_reduced_list.tst" );
gap> SetInfoLevel( InfoANUPQ, 1 );

# With automorphisms supplied, PqNextClass tested the exponent law only on
# a reduced list of words, which is enough only if the automorphisms induce
# all of GL(d, p) on a relatively free group; with fewer automorphisms the
# group came out too big (2^13 and 2^18 below instead of 2^4 and 2^12)
gap> procId := PqStart( FreeGroup( 4 ) : Prime := 2, Exponent := 2 );;
gap> PqPcPresentation( procId : ClassBound := 1 );;
gap> PqSupplyAutomorphisms( procId, [ IdentityMat( 4 ) ] );;
gap> PqNextClass( procId );;
gap> PqFactoredOrder( procId );
[ 2, 4 ]
gap> procId := PqStart( FreeGroup( 2 ) : Prime := 2, Exponent := 4 );;
gap> PqPcPresentation( procId : ClassBound := 1 );;
gap> PqSupplyAutomorphisms( procId, [ [ [ 0, 1 ], [ 1, 0 ] ] ] );;
gap> for class in [ 2 .. 6 ] do PqNextClass( procId ); od;
gap> PqFactoredOrder( procId );
[ 2, 12 ]

#
gap> PqQuitAll();
gap> STOP_TEST( "exponent_reduced_list.tst", 1 );
