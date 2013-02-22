% demonstrate_fetch.m
% Show how the data can be fetched from the Tunnel Data Server.

facility = 'T4'; shot_id = '7319'; channel = '110';

disp( 'First, get the metadata for the channel.' );
% and make use of the known order in which it is returned.  
[attrib, value] = fetch_channel_header( facility, shot_id, channel );
my_title = ['Facility ', facility, ', Shot ', char(value{9}), ...
            ', Channel ', char(value{14})];
my_xlabel = ['time in ', char(value{11})];
my_ylabel = [char(value{5}), ' in ', char(value{7})];

disp( 'Now, get the actual data and display graphically.' );
[t, v] = fetch_channel_data( facility, shot_id, channel );
subplot( 2, 1, 1 ); plot( t, v ); 
axis( [6000, 16000, 0, 40000] );
title( my_title ); xlabel( my_xlabel ); ylabel( my_ylabel );

disp( 'Filter and display again.' );
try
    % Attempt to use the signal-processing toolbox function.
    [b, a] = butter( 2, 0.05 );
catch
    % Set the filter coefficients manually.
    b = [0.0055427 0.0110854 0.0055427];
    a = [1.0 -1.7786318 0.8008026];
end
vf = filter( b, a, v );
subplot( 2, 1, 2 ); plot( t, vf ); 
axis( [6000, 16000, 0, 40000] );
xlabel( my_xlabel ); ylabel( my_ylabel );

disp( 'Averages during test-time.' );
tt = t > 7500 & t <= 8000;
m = sum( v .* tt ) / sum( tt );
mf = sum( vf .* tt ) / sum( tt );
disp( sprintf( 'Raw: %f, Filtered: %f', m, mf ) );
