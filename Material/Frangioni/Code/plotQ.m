function plotQ( Q , q , range )

f = @(x,y) 0.5 * [ x y ] * Q * [ x ; y ] + q' * [ x ; y ];

warning( 'off' , 'all' );

fcontour( f , range , 'LineColor' , 'k' , 'LineWidth' , 1 );

warning( 'on' , 'all' );

end