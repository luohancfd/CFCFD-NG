% test_data_server.m
% Show how the data can be fetched from the Tunnel Data Server.
facility = 'T4'; shot_id = '7319'; channel = '110';
useJava = 0; useCellArrays = 1;

% This script can be used in both MATLAB and Octave.
% Because there are differences, the script needs to be aware
% of its environment.
inOctave = exist( 'OCTAVE_VERSION' ) > 0;
if inOctave
    useJava = 0;       % The Java VM is only accessible from MATLAB 6.x.
    useCellArrays = 1; % list objects are deprecated in recent versions
end

% Here is the actual interaction with the Tunnel-Data Server.
tic;
if useJava
    disp( 'Talk to Tunnel Data Server via Java VM.' );
    disp( 'First, get the metadata.' );
    [attrib, value] = fetch_channel_header( facility, shot_id, channel );
    disp( 'Now, get the actual data.' );
    [t, v] = fetch_channel_data( facility, shot_id, channel );
else
    disp( 'Talk to the Tunnel Data Server via a system call.' );
    disp( 'First, get the metadata.' );
    [attrib, value] = get_data( facility, shot_id, channel, 'info' );
    disp( 'Now, get the actual data.' );
    [t, v] = get_data( facility, shot_id, channel, 'data' );
end
disp( sprintf( 'Elapsed time %f seconds.', toc ) );


% Build the plot titles and axis labels.
% Make use of the known order in which the metadata is returned.  
if useCellArrays
    % For MATLAB, metadata is in cell-arrays
    my_title = ['Facility ', facility, ', Shot ', char(value{9}), ...
                ', Channel ', char(value{14})];
    my_xlabel = ['time in ', char(value{11})];
    my_ylabel = [char(value{5}), ' in ', char(value{7})];
else
    % metadata is in list objects
    my_title = ['Facility ', facility, ', Shot ', nth(value,9), ...
                ', Channel ', nth(value,14)];
    my_xlabel = ['time in ', nth(value,11)];
    my_ylabel = [nth(value,5), ' in ', nth(value,7)];
endif

subplot( 2, 1, 1 ); 
if inOctave
    % set up annotation and then plot
    axis( [6000, 16000, 0, 40000] );
    title( [my_title, ', Raw Data'] ); 
    xlabel( my_xlabel ); ylabel( my_ylabel );
    plot( t, v );
else
    % first plot then annotate.
    plot( t, v );
    axis( [6000, 16000, 0, 40000] );
    title( [my_title, ', Raw Data'] ); 
    xlabel( my_xlabel ); ylabel( my_ylabel );
end


disp( 'Filter and plot data again.' );
try
    % Attempt to use the signal-processing toolbox function.
    [b, a] = butter( 2, 0.05 );
catch
    % But, if it is not available, set the filter coefficients manually.
    b = [0.0055427 0.0110854 0.0055427];
    a = [1.0 -1.7786318 0.8008026];
end
vf = filter( b, a, v );

subplot( 2, 1, 2 ); 
if inOctave
    % set up annotation and then plot
    axis( [6000, 16000, 0, 40000] );
    title( [my_title, ', Filtered Data'] ); 
    xlabel( my_xlabel ); ylabel( my_ylabel );
    plot( t, vf ); 
else
    % for MATLAB plot then annotate.
    plot( t, vf ); 
    axis( [6000, 16000, 0, 40000] );
    title( [my_title, ', Filtered Data'] ); 
    xlabel( my_xlabel ); ylabel( my_ylabel );
end

disp( 'Averages during test-time.' );
tt = t > 7500 & t <= 8000;
m = sum( v .* tt ) / sum( tt );
mf = sum( vf .* tt ) / sum( tt );
disp( sprintf( 'Raw: %f, Filtered: %f', m, mf ) );
