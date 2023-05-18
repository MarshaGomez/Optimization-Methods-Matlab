function fc( f , r )

warning( 'off' , 'all' );

fcontour( @(x,y) f( [ x ; y ] ) , r , 'LineColor' , 'k' , ...
                                      'LineWidth' , 1 );

ax = gca;
ax.FontSize = 16;
ax.Position = [ 0.05 0.05 0.94 0.94 ];
ax.Toolbar.Visible = 'off';

warning( 'on' , 'all' );

end

