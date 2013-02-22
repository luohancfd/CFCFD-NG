function content_text = fetch_text_from_server( ...
    facilityName, shotName, channelName, partName )
% FETCH_TEXT_FROM_SERVER fetches the requested data from the server.
% Input:
% facilityName : string, e.g. 'T3', 'T4', 'X1', 'X2' or 'X3'
% shotName     : string specifying the base-file-name for the shot.
% channelName  : string specifying the channel number.
% partName     : string, either 'data' or 'info', specifying
%                which part of the data file that we want to get.
% Output:
% Returns the data as a cell-arrays of strings.

% Peter J. 19-April-03

protocol = 'http';
% Set the following to appropriate values.
% host_name = '192.168.0.102';
host_name = 'www.mech.uq.edu.au';
with_passwd = 1;
user_name = 'tdguest';
passwd = 'tdpasswd';

% Set up the combined script-name and query string from the supplied pieces.
% It is not necessary to send shot and part names.
if isempty(partName) && isempty(channelName)
    cgi_script_name = sprintf( ...
        '/cgi-bin/tds/td_server.tcl?facility=%s+shot=%s', ...
        facilityName, shotName );
elseif isempty(partName)
    cgi_script_name = sprintf( ...
        '/cgi-bin/tds/td_server.tcl?facility=%s+shot=%s+channel=%s', ...
        facilityName, shotName, channelName );
else    
    cgi_script_name = sprintf( ...
        '/cgi-bin/tds/td_server.tcl?facility=%s+shot=%s+channel=%s+part=%s', ...
        facilityName, shotName, channelName, partName );
end

url = java.net.URL( protocol, host_name, cgi_script_name );
if with_passwd == 1
    % Encode the username and passwd before setting the Authorization
    connection = url.openConnection;
    text_string = [ user_name, ':', passwd ];
    % encoded_string = base64encode( text_string ); % not working
    encoded_string = 'dGRndWVzdDp0ZHBhc3N3ZA==';    % encoded manually
    connection.setRequestProperty( 'Authorization', ['Basic ', encoded_string] );
    connection.connect;
    is = connection.getInputStream;
else
    % We will operate without setting the Authorization property.     
    is = openStream( url );
end

isr = java.io.InputStreamReader( is );
ibr = java.io.BufferedReader( isr );

line_of_text = char( readLine( ibr ) );
count = 0;
content_text = {};  % Use a cell-array to collect the contents
while ~isempty(line_of_text)
    count = count + 1;
    content_text{count} = line_of_text;
    line_of_text = char( readLine( ibr ) );
end
