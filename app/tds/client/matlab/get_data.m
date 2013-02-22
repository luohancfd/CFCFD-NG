function [t, v] = get_data( facilityName, shotName, channelName, partName )
% GET_DATA gets specified data from the Tunnel-Data Server.
% Input:
% facilityName : string, e.g. 'T3', 'T4', 'X1', 'X2' or 'X3'
% shotName     : string specifying the base-file-name for the shot.
% channelName  : string specifying the channel number.
% partName     : string, either 'data' or 'info', specifying
%                which part of the data file that we want to get.
% Output:
% Depending on what is requested, the returned valuse may be:
% (1) the data as two numeric-arrays (vectors)
% (2) the metadata as two cell-arrays (or lists in Octave) of strings
% (3) the shot list as a cell-array (or lists in Octave) of strings
% (4) the channel names and Ids as two cell-arrays (or lists in Octave) of strings.

% PJ, 26-dec-01, 22-Apr-03

% First, decide what our action is to be for this call.
if strcmp( shotName, 'list' ) == 1
    action = 'getShotList';
    channelName = 'dummy';
    partName = 'dummy';
elseif strcmp( channelName, 'list' ) == 1
    action = 'getChannelList';
    partName = 'dummy';
elseif strcmp( partName, 'info' ) == 1
    action = 'getChannelInfo';
else
    action = 'getChannelData';
end

% select the mode of interaction with the HTTP/CGI server
%    0 = wget
%    1 = tcl_web_client 
useTclClientScript = 0;   
% Matlab won't know about the predefined constant OCTAVE_VERSION
inOctave = exist( 'OCTAVE_VERSION' ) > 0;
useCellArrays = 1;
if inOctave
    % Octave has the choice of cell arrays or lists, 
    % however, in recent versions of Octave (say >= 2.1.72)
    % list objects are deprecated.
    useCellArrays = 1;
endif

% Second, put together the command line options and then
% run an external process to get the data from the server.

if useTclClientScript == 1
    % If wget is not available, there is a TCL client script
    % that can talk to the Tunnel-Data Server via HTTP.
    cmd = ['./td_web_client.tcl ', facilityName, ' ', shotName, ' ', ...
           channelName, partName]
else
    % Set up a URL for wget to do the work. 
    server = 'www.mech.uq.edu.au';
    % server = 'galaxy5';  % test server
    userName = 'tdguest';
    password = 'tdpasswd';
    cgiCall = ['/cgi-bin/tds/td_server.tcl?facility=', facilityName, ...
               '+shot=', shotName, ...
               '+channel=', channelName, ...
               '+part=', partName];
    cmd = ['wget', ' --output-document=temporary.mat', ...
           ' --http-user=', userName, ' --http-passwd=', password, ...
           ' --quiet', ...
           ' http://', server, cgiCall];
end
% disp( cmd );
if inOctave
    % Use an Octave style system call.
    [result, status] = system(cmd);
else
    % Matlab style
    [status, result] = system(cmd);
end
if status ~= 0
    disp( 'The external process did not run successfully.' );
    disp( 'Maybe the server could not be found.' );
    t = [];
    v = [];
    return;
end 

% Before trying to pick up the requested data,
% peek at the first line of the temporary file 
% to see if the Tunnel-Data Server sent an error message.
[fid, msg] = fopen( 'temporary.mat', 'rt' );
line_of_text = fgetl( fid );
fclose( fid );
if ~isempty( findstr( line_of_text, 'td_server_error' ))
    disp( 'The data that was requested could not be found.' );
    disp( 'The server returned the following message:' );
    disp( line_of_text );
    t = [];
    v = [];
    return;
end

% Finally, pick up the local file containing the returned data
% and separate it into its components so that it can be returned.
if strcmp( action, 'getChannelData' ) == 1
    % The temporary file contains two columns of numbers.
    % Split into two numeric arrays (vectors) for return.
    if inOctave
        load temporary.mat;
    else
        load -ascii temporary.mat;
    end
    t = temporary(:,1);
    v = temporary(:,2);
elseif strcmp( action, 'getChannelInfo' ) == 1 || ...
       strcmp( action, 'getChannelList' ) == 1
    % The temporary file contains lines with two words each.
    count = 0; 
    if useCellArrays
        t = {}; v = {}; % use cell-arrays
    else
        t = list; v = list; % use lists
    end
    [fid, msg] = fopen( 'temporary.mat', 'rt' );
    line_of_text = fgetl( fid );
    while ~isnumeric( line_of_text )
        count = count + 1;
        if inOctave
            [a, b, n] = sscanf( line_of_text, '%s %s', 'C' );
        else
            [a, b] = strread( line_of_text, '%s %s' );
        end
	if useCellArrays
            t{count} = char(a);
            v{count} = char(b);
	else
            t = append( t, a );
            v = append( v, b );
	endif
        line_of_text = fgetl( fid );
    end
    fclose( fid );
elseif strcmp( action, 'getShotList' ) == 1
    % The temporary file contains a collection of shot names,
    % 10 per line, I believe.
    count = 0; 
    if useCellArrays
        t = {}; v = {}; % use cell-arrays
    else
        t = list; v = list; % use lists
    end
    [fid, msg] = fopen( 'temporary.mat', 'rt' );
    line_of_text = fgetl( fid );
    while ~isnumeric( line_of_text )
        count = count + 1;
        if inOctave
            t{count} = line_of_text;
        else
            t = append( t, line_of_text );
        end
        line_of_text = fgetl( fid );
    end
    fclose( fid );
end

